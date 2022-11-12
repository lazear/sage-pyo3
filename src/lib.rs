use annotate::AnnotatedPeak;
use psm::Psm;
use pyo3::exceptions::{PyFileNotFoundError, PyValueError};
use pyo3::prelude::*;
use sage_core::database::{Builder, IndexedDatabase, Parameters};
// use sage_core::fasta::Digest;
use sage_core::mass::{Mass, Residue, Tolerance};
use sage_core::peptide::Peptide;
use sage_core::scoring::Scorer;
use std::collections::HashMap;

mod annotate;
mod lfq;
mod psm;
mod spectra;

/// Python bindings to the Sage proteomic search engine
#[pymodule]
fn sage_proteomics(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<Database>()?;
    m.add_class::<psm::Psm>()?;
    m.add_class::<annotate::AnnotatedPeak>()?;
    let spectra_module = PyModule::new(_py, "spectra")?;
    spectra_module.add_class::<spectra::Mzml>()?;
    spectra_module.add_class::<spectra::Spectrum>()?;
    spectra_module.add_class::<spectra::Peak>()?;
    spectra_module.add_class::<spectra::Precursor>()?;
    m.add_submodule(spectra_module)?;
    Ok(())
}

#[pyclass]
struct Database {
    inner: IndexedDatabase,
    params: Parameters,
}

#[pymethods]
impl Database {
    /// new(fasta, /, decoy_prefix, static_mods, variable_mods)
    /// --
    /// Create a new Sage database
    #[new]
    // #[args(decoy_prefix="\"rev_\"", static_mods="{\"C\": 57.0215}")]
    #[args(decoy_tag = "\"rev_\".into()", generate_decoys = "true")]
    fn new(
        fasta: &str,
        decoy_tag: Option<String>,
        generate_decoys: Option<bool>,
        static_mods: Option<HashMap<char, f32>>,
        variable_mods: Option<HashMap<char, f32>>,
    ) -> PyResult<Self> {
        let mut builder = Builder {
            decoy_tag,
            generate_decoys,
            static_mods,
            variable_mods,
            ..Default::default()
        };
        builder.update_fasta(fasta.into());
        let params = builder.make_parameters();
        let db = params
            .clone()
            .build()
            .map_err(|e| PyErr::new::<PyFileNotFoundError, _>(fasta.to_string()))?;
        Ok(Self { inner: db, params })
    }

    #[getter]
    /// Number of fragment ions in database
    fn fragments(&mut self) -> usize {
        self.inner.fragments.len()
    }

    /// annotate_sequence(spectrum, sequence, / tolerance_ppm, charge, mods)
    /// --
    /// Annotate a MS2 spectrum with a provided peptide sequence
    fn annotate_sequence(
        &self,
        spectrum: spectra::Spectrum,
        sequence: String,
        tolerance_ppm: Option<f32>,
        charge: Option<u8>,
        mods: Option<HashMap<char, f32>>,
    ) -> PyResult<Vec<AnnotatedPeak>> {
        for c in sequence.chars() {
            if !sage_core::mass::VALID_AA.contains(&c) {
                return Err(PyErr::new::<PyValueError, _>(format!(
                    "Sequence {} contains invalid amino acid {}",
                    sequence, c
                )));
            }
        }
        let query = &spectrum.into();
        let sequence = sequence
            .chars()
            .map(|ch| Residue::Just(ch))
            .collect::<Vec<_>>();

        let mut peptide = Peptide {
            decoy: false,
            monoisotopic: sequence.iter().map(|m| m.monoisotopic()).sum::<f32>(),
            sequence,
            missed_cleavages: 0,
            nterm: None,
        };

        if let Some(mods) = mods {
            for (resi, mass) in mods {
                peptide.static_mod(resi, mass);
            }
        }

        Ok(annotate::annotate_peaks(
            query,
            peptide,
            tolerance_ppm.unwrap_or(10.0),
            charge,
        ))
    }

    /// annotate_psm(spectrum, psm, / tolerance_ppm, charge)
    /// --
    /// Annotate a MS2 spectrum with a provided PSM
    fn annotate_psm(
        &self,
        spectrum: spectra::Spectrum,
        psm: &Psm,
        tolerance_ppm: Option<f32>,
        charge: Option<u8>,
    ) -> Vec<AnnotatedPeak> {
        let query = &spectrum.into();
        let peptide = self.inner[psm.peptide_ix].clone();
        annotate::annotate_peaks(query, peptide, tolerance_ppm.unwrap_or(10.0), charge)
    }

    /// seach(spectrum, /, report_psms)
    /// --
    /// Search and score a single MS2 spectra, returning `report_psms` PSM objects
    #[args(report_psms = 1)]
    fn search(
        &self,
        spectrum: spectra::Spectrum,
        report_psms: Option<usize>,
    ) -> PyResult<Vec<Psm>> {
        if spectrum.level != 2 {
            return Err(PyErr::new::<PyValueError, _>(format!(
                "Scan {} has MS level = {}",
                spectrum.title, spectrum.level
            )));
        }
        let query = &spectrum.into();
        let scorer = Scorer::new(
            &self.inner,
            Tolerance::Ppm(-20.0, 20.0),
            Tolerance::Ppm(-10.0, 10.0),
            -1,
            3,
            Some(1),
            150.0,
            2000.0,
            false,
        );

        Ok(scorer
            .score(query, report_psms.unwrap_or(1))
            .into_iter()
            .map(Into::into)
            .collect())
    }

    /// seach(file, scan, /, report_psms)
    /// --
    /// Search and score a single MS2 spectra from a given file and scan ID,
    /// returning `report_psms` PSM objects
    fn search_by_id(
        &self,
        file: &spectra::Mzml,
        title: &str,
        report_psms: Option<usize>,
    ) -> PyResult<Vec<Psm>> {
        // let query = sage_core::spectrum::find_spectrum_by_id(&file.spectra, scan)
        // .ok_or(PyErr::new::<pyo3::exceptions::PyValueError, _>(scan))?;
        let query = file.get_sage_spectra(title)?;

        let scorer = Scorer::new(
            &self.inner,
            Tolerance::Ppm(-20.0, 20.0),
            Tolerance::Ppm(-10.0, 10.0),
            -1,
            3,
            Some(1),
            150.0,
            2000.0,
            false,
        );

        Ok(scorer
            .score(query, report_psms.unwrap_or(1))
            .into_iter()
            .map(Into::into)
            .collect())
    }
}
