// Copyright (c) 2020, David Michael Barr <b@rr-dav.id.au>. All rights reserved
//
// This source code is subject to the terms of the BSD 2 Clause License.

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
struct Vector {
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

#[derive(Copy, Clone)]
struct RotateZoom {
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

impl RotateZoom {
    fn new(width: i32, height: i32) -> Self {
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

fn main() {
    const WIDTH: usize = 1920;
    const HEIGHT: usize = 1080;
    let mut plane = Plane::<u8>::new(WIDTH, HEIGHT, 0, 0, 0, 0);
    let mut u = Vector { x: 1 << Q, y: 0 };
    let v = Vector { x: 16374, y: 257 };
    let rz = RotateZoom::new(WIDTH as i32, HEIGHT as i32);
    for frame in 0..1800 {
        fill_plane(&mut plane, rz * u, 8);
        u = u * v;
        image::GrayImage::from_fn(WIDTH as u32, HEIGHT as u32, |x, y| {
            image::Luma([plane.p(x as usize, y as usize)])
        })
        .save(format!("mandelbrot-{:03}.png", frame))
        .unwrap();
        println!("Frame {}", frame);
    }
}
