use pyo3::{exceptions::PyFileNotFoundError, prelude::*};
use rayon::prelude::*;
use sage_core::{
    mass::{Tolerance, PROTON},
    spectrum::{ProcessedSpectrum, SpectrumProcessor},
};
use std::collections::HashMap;

use crate::{
    lfq::{self, Xic},
    psm::Psm,
};

/// A Python module implemented in Rust.
#[pymodule]
fn spectra(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<Mzml>()?;
    m.add_class::<Spectrum>()?;
    m.add_class::<Peak>()?;
    m.add_class::<Precursor>()?;
    Ok(())
}

#[pyclass]
pub struct Ms2Iter {
    iter: Box<dyn Iterator<Item = Spectrum> + Send>,
}

#[pymethods]
impl Ms2Iter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<Spectrum> {
        slf.iter.next()
    }
}

#[pyclass]
pub struct Mzml {
    pub file: String,
    pub spectra: Vec<ProcessedSpectrum>,
    last_scan: usize,
    // Map spectrum title to index into `spectra` vector
    title_to_idx: HashMap<String, usize>,
}

impl Mzml {
    pub fn get_sage_spectra(&self, scan: &str) -> PyResult<&ProcessedSpectrum> {
        self.title_to_idx
            .get(scan)
            .and_then(|idx| self.spectra.get(*idx))
            .ok_or(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "invalid scan: {}",
                scan
            )))
    }
}

#[pymethods]
impl Mzml {
    #[new]
    pub fn new(path: &str) -> PyResult<Self> {
        let sp = SpectrumProcessor::new(150, 150.0, 2000.0, true, 0);
        let spectra = sage_cloudpath::read_mzml(path)
            .map_err(|e| PyErr::new::<PyFileNotFoundError, _>(path.to_string()))?
            .into_par_iter()
            .map(|spec| sp.process(spec))
            .collect::<Vec<_>>();
        let title_to_idx = spectra
            .iter()
            .enumerate()
            .map(|(idx, spec)| (spec.id.clone(), idx))
            .collect();

        Ok(Self {
            file: path.to_string(),
            spectra,
            last_scan: 0,
            title_to_idx,
        })
    }

    #[getter]
    pub fn scans(&self) -> usize {
        self.spectra.len()
    }

    pub fn __repr__(&self) -> PyResult<String> {
        Ok(format!("{} [{} scans]", self.file, self.spectra.len()))
    }

    pub fn get_spectrum(&self, scan: &str) -> PyResult<Spectrum> {
        self.get_sage_spectra(scan).cloned().map(Into::into)
    }

    /// Return MS1 peaks for extracted ion chromatogram
    #[args(rt_tolerance = "2.5")]
    pub fn xic(&self, psm: &Psm, charge: u8, rt_tolerance: Option<f32>) -> Vec<Xic> {
        let tol = rt_tolerance.unwrap_or(2.5);
        lfq::xic(
            &self.spectra,
            psm.rt - tol,
            psm.rt + tol,
            (psm.expmass - psm.isotope_error - PROTON) / charge as f32,
        )
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<Spectrum> {
        let x = slf.last_scan;
        slf.last_scan += 1;
        slf.spectra.get(x).cloned().map(Into::into)
    }
}

#[pyclass]
#[derive(Clone)]
pub struct Precursor {
    #[pyo3(get)]
    /// Precursor selected ion m/z
    pub mz: f32,
    #[pyo3(get)]
    /// Precursor selected ion intensity
    pub intensity: Option<f32>,
    #[pyo3(get)]
    /// Precursor selected ion charge
    pub charge: Option<u8>,
    #[pyo3(get)]
    /// Precursor scan number
    pub spectrum_ref: Option<String>,

    isolation_window: Option<Tolerance>,
}

#[pymethods]
impl Precursor {
    fn __repr__(&self) -> PyResult<String> {
        Ok(format!(
            "Precursor scan {:?}: m/z={}, z={:?}, int={:?}",
            self.spectrum_ref, self.mz, self.charge, self.intensity
        ))
    }
}

impl From<sage_core::spectrum::Precursor> for Precursor {
    fn from(val: sage_core::spectrum::Precursor) -> Self {
        Self {
            mz: val.mz,
            intensity: val.intensity,
            charge: val.charge,
            spectrum_ref: val.spectrum_ref,
            isolation_window: val.isolation_window,
        }
    }
}

impl Into<sage_core::spectrum::Precursor> for Precursor {
    fn into(self) -> sage_core::spectrum::Precursor {
        sage_core::spectrum::Precursor {
            mz: self.mz,
            intensity: self.intensity,
            charge: self.charge,
            spectrum_ref: self.spectrum_ref,
            isolation_window: self.isolation_window,
        }
    }
}

#[pyclass]
#[derive(Copy, Clone)]
pub struct Peak {
    #[pyo3(get)]
    /// Mass of peak (less proton), assumed to be z=1
    pub mass: f32,
    #[pyo3(get)]
    /// Peak intensity
    pub intensity: f32,
}

#[pymethods]
impl Peak {
    fn __repr(&self) -> String {
        format!(
            "MH+ = {}, Intensity = {}",
            self.mass + PROTON,
            self.intensity
        )
    }
}

impl From<sage_core::spectrum::Peak> for Peak {
    fn from(val: sage_core::spectrum::Peak) -> Self {
        Self {
            mass: val.mass,
            intensity: val.intensity,
        }
    }
}

impl Into<sage_core::spectrum::Peak> for Peak {
    fn into(self) -> sage_core::spectrum::Peak {
        sage_core::spectrum::Peak {
            mass: self.mass,
            intensity: self.intensity,
        }
    }
}

#[pyclass]
#[derive(Clone)]
pub struct Spectrum {
    /// MSn level
    #[pyo3(get)]
    pub level: u8,
    /// Scan ID
    #[pyo3(get)]
    pub title: String,
    /// Retention time
    #[pyo3(get)]
    pub scan_start_time: f32,
    /// Ion injection time
    #[pyo3(get)]
    pub ion_injection_time: f32,
    /// Selected ions for precursors, if `level > 1`
    #[pyo3(get)]
    pub precursors: Vec<Precursor>,
    /// MS peaks
    #[pyo3(get)]
    pub peaks: Vec<Peak>,
    #[pyo3(get)]
    /// Total MS2 intensity
    pub total_intensity: f32,
}

impl From<sage_core::spectrum::ProcessedSpectrum> for Spectrum {
    fn from(val: sage_core::spectrum::ProcessedSpectrum) -> Self {
        Self {
            level: val.level,
            title: val.id,
            scan_start_time: val.scan_start_time,
            ion_injection_time: val.ion_injection_time,
            precursors: val.precursors.into_iter().map(Into::into).collect(),
            peaks: val.peaks.into_iter().map(Into::into).collect(),
            total_intensity: val.total_intensity,
        }
    }
}

impl Into<sage_core::spectrum::ProcessedSpectrum> for Spectrum {
    fn into(self) -> sage_core::spectrum::ProcessedSpectrum {
        sage_core::spectrum::ProcessedSpectrum {
            level: self.level,
            id: self.title,
            file_id: 0,
            scan_start_time: self.scan_start_time,
            ion_injection_time: self.ion_injection_time,
            precursors: self.precursors.into_iter().map(Into::into).collect(),
            peaks: self.peaks.into_iter().map(Into::into).collect(),
            total_intensity: self.total_intensity,
        }
    }
}
