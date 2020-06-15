// Copyright (c) 2020, David Michael Barr <b@rr-dav.id.au>. All rights reserved
//
// This source code is subject to the terms of the BSD 2 Clause License.

use v_frame::frame::Frame;
use v_frame::plane::Plane;
use v_frame::prelude::*;

const Q: i32 = 14;
const P: i32 = 8;

fn mandelbrot_pixel(r: i32, i: i32) -> u16 {
    let mut z_r = 0i32;
    let mut z_i = 0i32;
    let mut z_rr = z_r * z_r;
    let mut z_ii = z_i * z_i;
    let low_bits = mix(r, i) & 63;
    for p in (1..=63).rev() {
        let z2_r = (z_rr - z_ii) >> Q;
        let z2_i = (z_r * z_i) >> (Q - 1);
        z_r = z2_r + r;
        z_i = z2_i + i;
        if let (Some(rr), Some(ii)) = (z_r.checked_mul(z_r), z_i.checked_mul(z_i)) {
            if let Some(norm) = rr.checked_add(ii) {
                if norm < (1 << (2 * Q + 2)) {
                    z_rr = rr;
                    z_ii = ii;
                    continue;
                }
            }
        }
        return (p << 6) | low_bits;
    }
    low_bits
}

fn mix(x: i32, y: i32) -> u16 {
    let h = x ^ (y >> 4) ^ (y << 4);
    (h ^ (h >> 7) ^ (h >> 4)) as u16
}

#[derive(Copy, Clone)]
pub struct Vector {
    x: i32,
    y: i32,
}

impl std::ops::Mul<Vector> for Vector {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        let ac = self.x.saturating_mul(rhs.x);
        let bd = self.y.saturating_mul(rhs.y);
        let a_b = self.x.saturating_add(self.y);
        let c_d = rhs.x.saturating_add(rhs.y);
        let abcd = a_b.saturating_mul(c_d);
        Self {
            x: ac.saturating_sub(bd) >> Q,
            y: abcd.saturating_sub(ac).saturating_sub(bd) >> Q,
        }
    }
}

impl Default for Vector {
    fn default() -> Self {
        Vector { x: 1 << Q, y: 0 }
    }
}

impl Vector {
    pub fn new(x: i32, y: i32, shift: usize) -> Self {
        assert!(shift <= Q as usize);
        Vector {
            x: x << (Q as usize - shift),
            y: y << (Q as usize - shift),
        }
    }

    pub fn pow(self, mut exponent: usize) -> Self {
        let mut product = if (exponent & 1) == 1 {
            self
        } else {
            Vector::default()
        };
        let mut factor = self;
        while exponent > 1 {
            exponent /= 2;
            factor = factor * factor;
            if (exponent & 1) == 1 {
                product = product * factor;
            }
        }
        product
    }
}

#[derive(Copy, Clone)]
pub struct RotateZoom {
    origin: Vector,
    col_step: Vector,
    row_step: Vector,
}

impl std::ops::Mul<Vector> for RotateZoom {
    type Output = Self;
    fn mul(self, rhs: Vector) -> Self {
        Self {
            origin: self.origin * rhs,
            col_step: self.col_step * rhs,
            row_step: self.row_step * rhs,
        }
    }
}

fn fill_slice<T: Pixel>(row: &mut [T], y: usize, rz: RotateZoom, bit_depth: usize) {
    for (x, v) in row.iter_mut().enumerate() {
        let r = ((x as i32 * rz.col_step.x + y as i32 * rz.row_step.x) >> P) + rz.origin.x;
        let i = ((x as i32 * rz.col_step.y + y as i32 * rz.row_step.y) >> P) + rz.origin.y;
        *v = T::cast_from(mandelbrot_pixel(r, i) >> (12 - bit_depth));
    }
}

fn fill_plane<T: Pixel>(plane: &mut Plane<T>, rz: RotateZoom, bit_depth: usize) {
    for (y, mut row) in plane
        .mut_slice(Default::default())
        .rows_iter_mut()
        .enumerate()
    {
        fill_slice(&mut row, y, rz, bit_depth);
    }
}

pub fn fill_frame<T: Pixel>(frame: &mut Frame<T>, rz: RotateZoom, bit_depth: usize) {
    for mut plane in frame.planes.iter_mut() {
        let xdec = plane.cfg.xdec;
        let ydec = plane.cfg.ydec;
        let rz = RotateZoom {
            origin: rz.origin,
            col_step: Vector {
                x: rz.col_step.x << xdec,
                y: rz.col_step.y << ydec,
            },
            row_step: Vector {
                x: rz.row_step.x << xdec,
                y: rz.row_step.y << ydec,
            },
        };
        fill_plane(&mut plane, rz, bit_depth);
    }
}

impl RotateZoom {
    pub fn new(width: usize, height: usize) -> Self {
        let width = width as i32;
        let height = height as i32;
        Self {
            origin: Vector {
                x: -(1 << (Q + 1)),
                y: -(height << Q) / width,
            },
            col_step: Vector {
                x: (1 << (P + Q + 1)) / width,
                y: 0,
            },
            row_step: Vector {
                x: 0,
                y: (1 << (P + Q + 1)) / width,
            },
        }
    }
}

#[allow(clippy::many_single_char_names)]
fn bt709_to_rgb<T: Pixel>(yuv: (T, T, T), bit_depth: usize) -> [T; 3] {
    const SHIFT: usize = 14;
    let y_den = 219 << (bit_depth - 8);
    let y_off = 16 << (bit_depth - 8);
    let uv_den = 224 << (bit_depth - 8);
    let uv_off = 128 << (bit_depth - 8);

    let y = (((i32::cast_from(yuv.0) - y_off) << SHIFT) + (y_den >> 1)) / y_den;
    let u = (((i32::cast_from(yuv.1) - uv_off) << SHIFT) + (uv_den >> 1)) / uv_den;
    let v = (((i32::cast_from(yuv.2) - uv_off) << SHIFT) + (uv_den >> 1)) / uv_den;

    let r = round_shift(y + ((20976 * v) >> SHIFT), SHIFT - bit_depth);
    let g = round_shift(y - ((3520 * u + 6236 * v) >> SHIFT), SHIFT - bit_depth);
    let b = round_shift(y + ((34865 * u) >> SHIFT), SHIFT - bit_depth);

    let r = r.min((1 << bit_depth) - 1).max(0);
    let g = g.min((1 << bit_depth) - 1).max(0);
    let b = b.min((1 << bit_depth) - 1).max(0);

    [T::cast_from(r), T::cast_from(g), T::cast_from(b)]
}

use std::path::PathBuf;
use structopt::StructOpt;

#[derive(Debug, StructOpt)]
struct Opts {
    #[structopt(short, long, default_value = "1920")]
    width: usize,
    #[structopt(short, long, default_value = "1080")]
    height: usize,
    #[structopt(short, long, default_value = ".", parse(from_os_str))]
    output_dir: PathBuf,
}

fn main() {
    let o = Opts::from_args();

    let mut frame = Frame::<u8>::new_with_padding(o.width, o.height, ChromaSampling::Cs420, 0);
    let v = Vector::new(16374, 257, 14);
    let rz = RotateZoom::new(o.width, o.height);

    for frame_number in 0..1800 {
        fill_frame(&mut frame, rz * v.pow(frame_number), 8);
        let xdec = frame.planes[1].cfg.xdec;
        let ydec = frame.planes[1].cfg.ydec;
        image::RgbImage::from_fn(o.width as u32, o.height as u32, |x, y| {
            image::Rgb(bt709_to_rgb(
                (
                    frame.planes[0].p(x as usize, y as usize),
                    frame.planes[1].p(x as usize >> xdec, y as usize >> ydec),
                    frame.planes[2].p(x as usize >> xdec, y as usize >> ydec),
                ),
                8,
            ))
        })
        .save(
            o.output_dir
                .join(format!("mandelbrot-{:03}.png", frame_number)),
        )
        .unwrap();
        println!("Frame {}", frame_number);
    }
}
