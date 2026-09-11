//! Compare a completed output with its input on identical binned axes and scale.
use eframe::egui::{self, Color32, ColorImage, TextureHandle};
use std::path::PathBuf;
use std::sync::mpsc::{Receiver, channel};
const WIDTH: usize = 256;
const HEIGHT: usize = 128;
type Points = Vec<(u32, u32, u32)>;
struct Loaded {
    images: Vec<ColorImage>,
    description: String,
}
pub struct Preview {
    input: PathBuf,
    output: PathBuf,
    frame: usize,
    open: bool,
    rx: Option<Receiver<Result<Loaded, String>>>,
    textures: Vec<TextureHandle>,
    description: String,
}
impl Preview {
    pub fn new(input: PathBuf, output: PathBuf) -> Self {
        let mut preview = Self {
            input,
            output,
            frame: 1,
            open: true,
            rx: None,
            textures: vec![],
            description: String::new(),
        };
        preview.load();
        preview
    }
    fn load(&mut self) {
        let (tx, rx) = channel();
        self.rx = Some(rx);
        self.textures.clear();
        let input = self.input.clone();
        let output = self.output.clone();
        let index = self.frame.saturating_sub(1);
        std::thread::spawn(move || {
            let result = (|| {
                let (scans, raw) =
                    dnoise::validation::read_frame(&input, index).map_err(|e| e.to_string())?;
                let (out_scans, processed) =
                    dnoise::validation::read_frame(&output, index).map_err(|e| e.to_string())?;
                if scans != out_scans {
                    return Err("input/output scan counts differ".into());
                }
                Ok(heatmaps(scans, &raw, &processed))
            })();
            let _ = tx.send(result);
        });
    }
    pub fn show(&mut self, ctx: &egui::Context) {
        if let Some(rx) = &self.rx
            && let Ok(result) = rx.try_recv()
        {
            self.rx = None;
            match result {
                Ok(loaded) => {
                    self.description = loaded.description;
                    self.textures = loaded
                        .images
                        .into_iter()
                        .enumerate()
                        .map(|(i, image)| {
                            ctx.load_texture(
                                format!("comparison-{i}"),
                                image,
                                egui::TextureOptions::NEAREST,
                            )
                        })
                        .collect();
                }
                Err(e) => self.description = e,
            }
        }
        let mut open = self.open;
        egui::Window::new("Frame comparison").open(&mut open).default_width(920.).show(ctx,|ui| {
            ui.horizontal(|ui| {
                ui.label("Frame (1-based):");ui.add(egui::DragValue::new(&mut self.frame).range(1..=usize::MAX));
                if ui.add_enabled(self.rx.is_none(),egui::Button::new("Load frame")).clicked() {self.load();}
            });
            if self.rx.is_some() {ui.spinner();ctx.request_repaint_after(std::time::Duration::from_millis(100));}
            let width=(ui.available_width()-24.)/3.;
            ui.horizontal(|ui| {
                for (label,texture) in ["Before","After","Removed intensity"].into_iter().zip(&self.textures) {
                    ui.vertical(|ui| {ui.strong(label);ui.add(egui::Image::new(texture).fit_to_exact_size(egui::vec2(width,220.)));});
                }
            });
            ui.label(&self.description);
            ui.weak("Same axes and logarithmic intensity scale. Pixels sum points; removed intensity is the positive before-minus-after difference. This view does not establish scientific fidelity.");
        });
        self.open = open;
    }
}
fn heatmaps(scans: usize, raw: &Points, processed: &Points) -> Loaded {
    let max_tof = raw
        .iter()
        .chain(processed)
        .map(|p| p.1)
        .max()
        .unwrap_or(1)
        .max(1);
    let min_tof = raw.iter().chain(processed).map(|p| p.1).min().unwrap_or(0);
    let span = (max_tof - min_tof).max(1) as f64;
    let bin = |points: &Points| {
        let mut bins = vec![0f64; WIDTH * HEIGHT];
        for &(scan, tof, intensity) in points {
            let x = (((tof - min_tof) as f64 / span) * (WIDTH - 1) as f64).round() as usize;
            let y = ((scan as f64 / scans.saturating_sub(1).max(1) as f64) * (HEIGHT - 1) as f64)
                .round() as usize;
            if x < WIDTH && y < HEIGHT {
                bins[(HEIGHT - 1 - y) * WIDTH + x] += intensity as f64;
            }
        }
        bins
    };
    let before = bin(raw);
    let after = bin(processed);
    let removed: Vec<_> = before
        .iter()
        .zip(&after)
        .map(|(a, b)| (a - b).max(0.))
        .collect();
    let maximum = before
        .iter()
        .chain(&after)
        .copied()
        .fold(1., f64::max)
        .ln_1p();
    let images = [before, after, removed]
        .into_iter()
        .map(|bins| {
            let mut image = ColorImage::new([WIDTH, HEIGHT], vec![Color32::BLACK; WIDTH * HEIGHT]);
            for (pixel, value) in image.pixels.iter_mut().zip(bins) {
                let t = (value.ln_1p() / maximum).clamp(0., 1.);
                *pixel =
                    Color32::from_rgb((255. * t) as u8, (220. * t.sqrt()) as u8, (100. * t) as u8);
            }
            image
        })
        .collect();
    let total: f64 = raw.iter().map(|p| p.2 as f64).sum();
    let kept: f64 = processed.iter().map(|p| p.2 as f64).sum();
    Loaded {
        images,
        description: format!(
            "TOF index {min_tof}–{max_tof} (left → right); mobility scan 0–{} (bottom → top). {} → {} points; {:.1}% intensity retained.",
            scans.saturating_sub(1),
            raw.len(),
            processed.len(),
            if total > 0. { 100. * kept / total } else { 0. }
        ),
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn unchanged_frame_has_black_removed_panel_and_equal_panels() {
        let points = vec![(1, 10, 42), (2, 15, 100)];
        let loaded = heatmaps(3, &points, &points);
        assert_eq!(loaded.images[0].pixels, loaded.images[1].pixels);
        assert!(loaded.images[2].pixels.iter().all(|p| *p == Color32::BLACK));
    }
}
