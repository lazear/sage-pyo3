use pyo3::prelude::*;
use rayon::prelude::*;
use sage_core::database::binary_search_slice;
use sage_core::ion_series::{IonSeries, Kind};
use sage_core::mass::Tolerance;
use sage_core::peptide::Peptide;
use sage_core::spectrum::ProcessedSpectrum;

#[pyclass]
pub struct AnnotatedPeak {
    #[pyo3(get)]
    mass: f32,
    #[pyo3(get)]
    intensity: f32,
    #[pyo3(get)]
    ion: char,
    #[pyo3(get)]
    index: usize,
    #[pyo3(get)]
    charge: u8,
}

#[pymethods]
impl AnnotatedPeak {
    fn __repr__(&self) -> String {
        format!("{}{} {}+", self.ion, self.index, self.charge)
    }
}

/// Calculate full hyperscore for a given PSM
pub fn annotate_peaks(
    query: &ProcessedSpectrum,
    peptide: Peptide,
    tolerance_ppm: f32,
    charge: Option<u8>,
) -> Vec<AnnotatedPeak> {
    // Regenerate theoretical ions
    let mut fragments = IonSeries::new(&peptide, Kind::B)
        .enumerate()
        .map(|(idx, ion)| (idx + 1, ion))
        .chain(
            IonSeries::new(&peptide, Kind::Y)
                .enumerate()
                .map(|(idx, ion)| (peptide.sequence.len().saturating_sub(1 + idx), ion)),
        )
        .collect::<Vec<_>>();

    fragments.sort_unstable_by(|a, b| a.1.monoisotopic_mass.total_cmp(&b.1.monoisotopic_mass));

    let charge = charge.unwrap_or(2);
    query
        .peaks
        .par_iter()
        .flat_map_iter(|peak| {
            (1..charge).flat_map(|charge| {
                let peak = *peak;
                let mass = peak.mass * charge as f32;
                let (lo, hi) = Tolerance::Ppm(-tolerance_ppm, tolerance_ppm).bounds(mass);
                let window = binary_search_slice(
                    &fragments,
                    |frag, mz| frag.1.monoisotopic_mass.total_cmp(mz),
                    lo,
                    hi,
                );

                // for frag in fragments[window.0..window.1]
                (&fragments[window.0..window.1])
                    .iter()
                    .filter(move |frag| {
                        frag.1.monoisotopic_mass >= lo && frag.1.monoisotopic_mass <= hi
                    })
                    .map(move |(index, frag)| {
                        let ion = match frag.kind {
                            Kind::B => 'b',
                            Kind::Y => 'y',
                        };

                        AnnotatedPeak {
                            mass: peak.mass,
                            intensity: peak.intensity,
                            ion,
                            index: *index,
                            charge,
                        }
                    })
            })
        })
        .collect()
}
