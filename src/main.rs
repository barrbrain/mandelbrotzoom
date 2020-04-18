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

fn fill_slice(p: &mut [u8], stride: usize) {
    let y_step = (stride << (Q + 1)) / p.len();
    let x_step = (1 << (Q + 1)) / stride;
    for (y, row) in p.chunks_mut(stride).enumerate() {
        for (x, v) in row.iter_mut().enumerate() {
            let r = (x * x_step) as i32 - (1 << (Q + 1));
            let i = (y * y_step) as i32 - (1 << Q);
            *v = mandelbrot_pixel(r, i);
        }
    }
}

fn main() {
    const WIDTH: usize = 1920;
    const HEIGHT: usize = 1080;
    let mut buf = [0; WIDTH * HEIGHT];
    fill_slice(&mut buf, WIDTH);
    image::GrayImage::from_fn(WIDTH as u32, HEIGHT as u32, |x, y| {
        image::Luma([buf[(y * (WIDTH as u32) + x) as usize]])
    })
    .save("mandelbrot.png")
    .unwrap();
}
