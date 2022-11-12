use pyo3::prelude::*;
use sage_core::database::PeptideIx;
use sage_core::scoring::Feature;

#[pyclass]
pub struct Psm {
    pub peptide_ix: PeptideIx,
    /// Peptide sequence, including modifications e.g.: NC(+57.021)HK
    #[pyo3(get)]
    pub peptide: String,
    /// Peptide length
    #[pyo3(get)]
    pub peptide_len: usize,
    /// Proteins containing this peptide sequence
    #[pyo3(get)]
    pub proteins: Vec<String>,
    /// Number of proteins assigned to this peptide sequence
    #[pyo3(get)]
    pub num_proteins: usize,
    /// Spectrum title
    #[pyo3(get)]
    pub spectrum_title: String,
    /// Target/Decoy label, -1 is decoy, 1 is target
    #[pyo3(get)]
    pub decoy: bool,
    /// Experimental mass MH+
    #[pyo3(get)]
    pub expmass: f32,
    /// Calculated mass, MH+
    #[pyo3(get)]
    pub calcmass: f32,
    /// Reported precursor charge
    #[pyo3(get)]
    pub charge: u8,
    /// Retention time
    #[pyo3(get)]
    pub rt: f32,
    /// Difference between expmass and calcmass
    #[pyo3(get)]
    pub delta_mass: f32,
    /// C13 isotope error
    #[pyo3(get)]
    pub isotope_error: f32,
    /// Average ppm delta mass for matched fragments
    #[pyo3(get)]
    pub average_ppm: f32,
    /// X!Tandem hyperscore
    #[pyo3(get)]
    pub hyperscore: f64,
    /// Difference between hyperscore of this candidate, and the next best candidate
    #[pyo3(get)]
    pub delta_hyperscore: f64,
    /// Number of matched theoretical fragment ions
    #[pyo3(get)]
    pub matched_peaks: u32,
    /// Longest b-ion series
    #[pyo3(get)]
    pub longest_b: u32,
    /// Longest y-ion series
    #[pyo3(get)]
    pub longest_y: u32,
    #[pyo3(get)]
    pub longest_y_pct: f32,
    /// Number of missed cleavages
    #[pyo3(get)]
    pub missed_cleavages: u8,
    /// Fraction of matched MS2 intensity
    #[pyo3(get)]
    pub matched_intensity_pct: f32,
    /// Number of scored candidates for this spectrum
    #[pyo3(get)]
    pub scored_candidates: u32,
    /// Probability of matching exactly N peaks across all candidates Pr(x=k)
    #[pyo3(get)]
    pub poisson: f64,
}

#[pymethods]
impl Psm {
    fn __repr__(&self) -> String {
        format!("{} {:?}", self.peptide, self.proteins)
    }
}

impl From<Feature> for Psm {
    fn from(p: Feature) -> Self {
        Self {
            peptide_ix: p.peptide_idx,
            peptide: p.peptide,
            peptide_len: p.peptide_len,
            proteins: p.proteins.split(";").map(Into::into).collect(),
            num_proteins: p.num_proteins,
            spectrum_title: p.spec_id,
            decoy: p.label == -1,
            expmass: p.expmass,
            calcmass: p.calcmass,
            charge: p.charge,
            rt: p.rt,
            delta_mass: p.delta_mass,
            isotope_error: p.isotope_error,
            average_ppm: p.average_ppm,
            hyperscore: p.hyperscore,
            delta_hyperscore: p.delta_hyperscore,
            matched_peaks: p.matched_peaks,
            longest_b: p.longest_b,
            longest_y: p.longest_y,
            longest_y_pct: p.longest_y_pct,
            missed_cleavages: p.missed_cleavages,
            matched_intensity_pct: p.matched_intensity_pct,
            scored_candidates: p.scored_candidates,
            poisson: p.poisson,
        }
    }
}
