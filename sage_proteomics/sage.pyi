from typing import Optional, Dict, List
from . import spectra

class AnnotatedPeak:
    mass: float
    """Sage internally uses masses instead of m/z - peaks are assumed to be z=1"""
    intensity: float
    ion: str
    """b or y"""
    index: int
    """b1, b2, y3, y4, etc"""
    charge: int

class Psm:
    peptide: str
    """Peptide sequence, in ProForma notation"""
    peptide_len: int
    """Peptide length"""
    proteins: List[str]
    """Proteins from which this peptide may have originated"""
    num_proteins: int
    spectrum_title: str
    """MS2 spectrum title"""
    decoy: bool
    """Is this PSM a decoy?"""
    expmass: float
    """Experimental mass"""
    calcmass: float
    """Calculated monoisotopic mass"""
    charge: int
    """Observed or calculated precursor charge state"""
    rt: float
    """Retention time"""
    delta_mass: float
    """Difference between `expmass` and `calcmass` in ppm"""
    isotope_error: float
    """Isotopic error"""
    average_ppm: float
    """Average difference between experimental and theoretical MS2 peaks, in ppm"""
    hyperscore: float
    """X!Tandem hyperscore"""
    matched_peaks: int
    """Number of matched peaks"""
    longest_b: int
    """Longest continuous ladder of b-ion fragments"""
    longest_y: int
    """Longest continuous ladder of y-ion fragments"""
    longest_y_pct: int
    """Longest continuous ladder of y-ion fragments, as a % of peptide length"""
    missed_cleavages: int
    """Number of missed cleavages"""
    matched_intensity_pct: float
    """Percentage of explained MS2 intensity"""
    scored_candidates: int
    """Number of candidate peptides scored for this spectrum"""
    poisson: float
    """log10 probability of matching this many peaks across all candidates"""

class Database:
    """
    A class representating a FASTA database that has been digested and indexed
    """

    def __init__(
        self,
        path: str,
        decoy_tag: Optional[str] = "rev_",
        generate_decoys: Optional[bool] = True,
        static_mods: Optional[Dict[str, float]] = None,
        variable_mods: Optional[Dict[str, float]] = None,
    ) -> None:
        """
        Create a new Sage database
        """
    def annotate_sequence(
        self,
        spectrum: spectra.Spectrum,
        sequence: str,
        tolerance_ppm: Optional[float],
        charge: Optional[int],
        mods: Optional[Dict[str, float]],
    ) -> List[AnnotatedPeak]:
        """
        Annotated a MS2 `Spectrum`, returning a list of matching MS2 peaks
        """
    def annotate_psm(
        self,
        spectrum: spectra.Spectrum,
        psm: Psm,
        tolerance_ppm: Optional[float],
        charge: Optional[int],
    ) -> List[AnnotatedPeak]:
        """
        Given a peptide-spectrum match, return a list of the matched peaks
        """
    def search(
        self, spectrum: spectra.Spectrum, report_psms: Optional[int] = 1
    ) -> List[Psm]:
        """
        Search and score a single MS2 spectra, returning a list containing
        `report_psms` Psms
        """
    def search_by_id(
        self, file, title: str, report_psms: Optional[int] = 1
    ) -> List[Psm]:
        """
        Locate and score a single MS2 spectra from `file` identified by spectrum title.
        Return a list of Psms of length `report_psms`
        """
