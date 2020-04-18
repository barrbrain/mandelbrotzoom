// Copyright (c) 2020, David Michael Barr <b@rr-dav.id.au>. All rights reserved
//
// This source code is subject to the terms of the BSD 2 Clause License.

const Q: i32 = 14;

fn mandelbrot_pixel(r: i32, i: i32) -> u8 {
    let mut z_r = 0i32;
    let mut z_i = 0i32;
    let mut z_rr = z_r * z_r;
    let mut z_ii = z_i * z_i;
    for p in (1..=255).rev() {
        let z2_r = (z_rr - z_ii) >> Q;
        let z2_i = (z_r * z_i) >> (Q - 1);
        z_r = z2_r + r;
        z_i = z2_i + i;
        z_rr = z_r.saturating_mul(z_r);
        z_ii = z_i.saturating_mul(z_i);
        if z_rr.saturating_add(z_ii) >= (1 << (2 * Q + 2)) {
            return p;
        }
    }
    0
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

fn fill_slice(p: &mut [u8], stride: usize, rz: RotateZoom) {
    for (y, row) in p.chunks_mut(stride).enumerate() {
        for (x, v) in row.iter_mut().enumerate() {
            let r = x as i32 * rz.col_step.x + y as i32 * rz.row_step.x + rz.origin.x;
            let i = x as i32 * rz.col_step.y + y as i32 * rz.row_step.y + rz.origin.y;
            *v = mandelbrot_pixel(r, i);
        }
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
                x: (1 << (Q + 1)) / width,
                y: 0,
            },
            row_step: Vector {
                x: 0,
                y: (1 << (Q + 1)) / width,
            },
        }
    }
}

fn main() {
    const WIDTH: usize = 1920;
    const HEIGHT: usize = 1080;
    let mut buf = [0; WIDTH * HEIGHT];
    let v = Vector {
        x: 254 << (Q - 8),
        y: 3 << (Q - 8),
    };
    let rz = RotateZoom::new(WIDTH as i32, HEIGHT as i32) * v;
    fill_slice(&mut buf, WIDTH, rz);
    image::GrayImage::from_fn(WIDTH as u32, HEIGHT as u32, |x, y| {
        image::Luma([buf[(y * (WIDTH as u32) + x) as usize]])
    })
    .save("mandelbrot.png")
    .unwrap();
}
