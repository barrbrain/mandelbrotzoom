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
    plane
        .mut_slice(Default::default())
        .rows_iter_mut()
        .enumerate()
        .collect::<Vec<_>>()
        .into_par_iter()
        .for_each(|(y, mut row)| {
            fill_slice(&mut row, y, rz, bit_depth);
        });
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

fn parse_chroma_sampling(s: &str) -> Result<ChromaSampling, String> {
    match s {
        "4:2:0" => Ok(ChromaSampling::Cs420),
        "4:2:2" => Ok(ChromaSampling::Cs422),
        "4:4:4" => Ok(ChromaSampling::Cs444),
        other => Err(format!("Unsupported Chroma Sampling {}", other)),
    }
}

fn parse_bit_depth(s: &str) -> Result<usize, String> {
    use std::str::FromStr;

    let v = usize::from_str(s).map_err(|e| e.to_string())?;

    if v != 8 && v != 10 && v != 12 {
        Err(format!("Unsupported bit_depth {}", v))
    } else {
        Ok(v)
    }
}

#[derive(Debug, StructOpt)]
struct Opts {
    #[structopt(short, long, default_value = "1920")]
    width: usize,
    #[structopt(short, long, default_value = "1080")]
    height: usize,
    #[structopt(short, long, default_value = "mandelbrot.y4m", parse(from_os_str))]
    output: PathBuf,
    #[structopt(short, long, default_value = "1800")]
    frames: usize,
    #[structopt(short, long, default_value = "4:2:0", parse(try_from_str = parse_chroma_sampling))]
    chroma_sampling: ChromaSampling,
    // #[structopt(short, long, default_value = "8", parse(try_from_str = parse_bit_depth))]
    // bit_depth: usize,
    #[structopt(long, default_value = "60")]
    fps: usize,
}

use rayon::prelude::*;
use std::io::Write;

fn write_y4m_frame<T: Pixel, W: Write>(y4m_enc: &mut y4m::Encoder<W>, rec: &Frame<T>, o: &Opts) {
    use std::slice;

    let bit_depth = 8;

    let planes = if o.chroma_sampling == ChromaSampling::Cs400 {
        1
    } else {
        3
    };
    let bytes_per_sample = if bit_depth > 8 { 2 } else { 1 };
    let (chroma_width, chroma_height) = o.chroma_sampling.get_chroma_dimensions(o.width, o.height);
    let pitch_y = o.width * bytes_per_sample;
    let pitch_uv = chroma_width * bytes_per_sample;

    let (mut rec_y, mut rec_u, mut rec_v) = (
        vec![128u8; pitch_y * o.height],
        vec![128u8; pitch_uv * chroma_height],
        vec![128u8; pitch_uv * chroma_height],
    );

    let (stride_y, stride_u, stride_v) = (
        rec.planes[0].cfg.stride,
        rec.planes[1].cfg.stride,
        rec.planes[2].cfg.stride,
    );

    for (line, line_out) in rec.planes[0]
        .data_origin()
        .chunks(stride_y)
        .zip(rec_y.chunks_mut(pitch_y))
    {
        if bit_depth > 8 {
            unsafe {
                line_out.copy_from_slice(slice::from_raw_parts::<u8>(
                    line.as_ptr() as *const u8,
                    pitch_y,
                ));
            }
        } else {
            line_out.copy_from_slice(
                &line.iter().map(|&v| u8::cast_from(v)).collect::<Vec<u8>>()[..pitch_y],
            );
        }
    }

    if planes > 1 {
        for (line, line_out) in rec.planes[1]
            .data_origin()
            .chunks(stride_u)
            .zip(rec_u.chunks_mut(pitch_uv))
        {
            if bit_depth > 8 {
                unsafe {
                    line_out.copy_from_slice(slice::from_raw_parts::<u8>(
                        line.as_ptr() as *const u8,
                        pitch_uv,
                    ));
                }
            } else {
                line_out.copy_from_slice(
                    &line.iter().map(|&v| u8::cast_from(v)).collect::<Vec<u8>>()[..pitch_uv],
                );
            }
        }
        for (line, line_out) in rec.planes[2]
            .data_origin()
            .chunks(stride_v)
            .zip(rec_v.chunks_mut(pitch_uv))
        {
            if bit_depth > 8 {
                unsafe {
                    line_out.copy_from_slice(slice::from_raw_parts::<u8>(
                        line.as_ptr() as *const u8,
                        pitch_uv,
                    ));
                }
            } else {
                line_out.copy_from_slice(
                    &line.iter().map(|&v| u8::cast_from(v)).collect::<Vec<u8>>()[..pitch_uv],
                );
            }
        }
    }

    let rec_frame = y4m::Frame::new([&rec_y, &rec_u, &rec_v], None);
    y4m_enc.write_frame(&rec_frame).unwrap();
}

use anyhow::Result;

fn main() -> Result<()> {
    let o = Opts::from_args();

    let v = Vector::new(16374, 257, 14);
    let rz = RotateZoom::new(o.width, o.height);

    let bit_depth = 8;

    let csp = match (o.chroma_sampling, bit_depth) {
        (ChromaSampling::Cs420, 8) => y4m::Colorspace::C420,
        (ChromaSampling::Cs420, 10) => y4m::Colorspace::C420p10,
        (ChromaSampling::Cs420, 12) => y4m::Colorspace::C420p12,
        (ChromaSampling::Cs422, 8) => y4m::Colorspace::C422,
        (ChromaSampling::Cs422, 10) => y4m::Colorspace::C422p10,
        (ChromaSampling::Cs422, 12) => y4m::Colorspace::C422p12,
        (ChromaSampling::Cs444, 8) => y4m::Colorspace::C444,
        (ChromaSampling::Cs444, 10) => y4m::Colorspace::C444p10,
        (ChromaSampling::Cs444, 12) => y4m::Colorspace::C444p12,
        _ => unreachable!(),
    };

    let out = std::fs::File::create(&o.output)?;

    let mut y4m_enc = y4m::encode(o.width, o.height, y4m::Ratio::new(o.fps, 1))
        .with_colorspace(csp)
        .write_header(out)
        .unwrap();
    let mut frame = Frame::<u8>::new_with_padding(o.width, o.height, ChromaSampling::Cs420, 0);

    (0..o.frames).for_each(|frame_number| {
        fill_frame(&mut frame, rz * v.pow(frame_number), 8);

        write_y4m_frame(&mut y4m_enc, &frame, &o);
    });

    Ok(())
}
