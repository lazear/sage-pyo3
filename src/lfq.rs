use crate::spectra::Peak;
use pyo3::prelude::*;
use rayon::prelude::*;
use sage_core::{database::binary_search_slice, mass::Tolerance, spectrum::ProcessedSpectrum};

#[pyclass]
pub struct Xic {
    #[pyo3(get)]
    rt: f32,
    #[pyo3(get)]
    mass: f32,
    #[pyo3(get)]
    intensity: f32,
}

pub fn xic(spectra: &[ProcessedSpectrum], rt_min: f32, rt_max: f32, mass: f32) -> Vec<Xic> {
    let (lo, hi) = Tolerance::Ppm(-5.0, 5.0).bounds(mass);

    let window = binary_search_slice(
        &spectra,
        |spec, rt| spec.scan_start_time.total_cmp(rt),
        rt_min,
        rt_max,
    );
    spectra[window.0..window.1]
        .par_iter()
        .filter(|s| s.level == 1)
        .flat_map_iter(|spectrum| {
            spectrum
                .peaks
                .iter()
                .filter(|peak| peak.mass >= lo && peak.mass <= hi)
                .map(|peak| Xic {
                    rt: spectrum.scan_start_time,
                    mass: peak.mass,
                    intensity: peak.intensity,
                })
        })
        .collect()
}
