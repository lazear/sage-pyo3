from typing import Optional, Dict, List

class Peak:
    mass: float
    """Sage internally uses masses instead of m/z - peaks are assumed to be z=1"""
    intensity: float

class Precursor:
    mz: float
    intensity: Optional[float]
    charge: Optional[int]
    spectrum_ref: Optional[str]

class Spectrum:
    level: int
    """MS level (e.g. 1, 2, 3)"""
    title: str
    """Spectrum title"""
    scan_start_time: float
    ion_injection_time: float
    precursors: List[Precursor]
    """Precursor information"""
    peaks: List[Peak]
    total_intensity: float

class Mzml:
    def __init__(self, path: str) -> None:
        """
        Read an mzML file

        Sage supports reading mzML files that are packaged in gzip containers,
        or files directly from AWS S3
        """
    def get_spectrum(self, title: str) -> Spectrum:
        """
        Find a spectrum by title
        """
    scans: int
    """Number of spectra in this file"""
