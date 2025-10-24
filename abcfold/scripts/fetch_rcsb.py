"""Module for downloading structure files from RCSB."""

import logging
from pathlib import Path

import requests
from tqdm import tqdm

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# Adapted from antid.io.struct.RCSBDownloader
class RCSBDownloader:
    """Download structure files from RCSB."""

    def __init__(
        self,
        out_dir: str | Path,
        make_subdir: bool = False,
        req_session: requests.Session | None = None,
        timeout: int = 10,
    ):
        """Initialize the downloader.

        Args:
            out_dir: Directory to save downloaded files.
            make_subdir: If True, always use the two characters in the middle of the PDB
                ID as subdirectory names. This is useful when downloading a large number
                of PDB files to avoid too many files in a single directory.
            req_session: Optional requests.Session object.
            timeout: Timeout for requests in seconds.

        """
        self.out_dir = Path(out_dir).expanduser().resolve()
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.subdir = make_subdir
        self.session = requests.Session() if req_session is None else req_session
        self.timeout = timeout

    def fetch_pdb(
        self, pdb_id: str, file_type: str = "bio", fallback_to_cif: bool = True
    ) -> Path:
        """Download a PDB file from RCSB.

        Args:
            pdb_id: The PDB ID.
            file_type: Bio-assembly1 (bio) or asymmetric unit (asu).
            fallback_to_cif: If True, fall back to downloading the mmCIF file if the PDB
            file is not found. For details, see https://www.rcsb.org/docs/general-help/structures-without-legacy-pdb-format-files

        Returns:
            Path to the downloaded file.

        """
        pdb_id = pdb_id.upper()
        file_type = self._check_file_type(file_type)
        out_dir = (
            Path(self.out_dir / file_type / pdb_id[-3:-1])
            if self.subdir
            else self.out_dir / file_type
        )
        out_dir.mkdir(parents=True, exist_ok=True)
        gz_pdb_file = out_dir / f"{pdb_id}.pdb.gz"
        if gz_pdb_file.exists():
            return gz_pdb_file

        pdb_url = (
            f"https://files.rcsb.org/download/{pdb_id}.pdb1.gz"
            if file_type == "bio"
            else f"https://files.rcsb.org/download/{pdb_id}.pdb.gz"
        )
        r = self.session.get(pdb_url, timeout=self.timeout)
        if r.status_code == 404 and fallback_to_cif:
            # Try mmCIF file if the PDB is nonexistent
            logger.warning(f"PDB for {pdb_id} not found, falling back to mmCIF.")
            return self.fetch_mmcif(pdb_id)

        r.raise_for_status()
        with open(gz_pdb_file, "wb") as f:
            f.write(r.content)
        return gz_pdb_file

    def fetch_mmcif(self, pdb_id: str, file_type: str = "bio") -> Path:
        """Download a mmCIF file from RCSB.

        Args:
            pdb_id: The PDB ID.
            file_type: Bio-assembly1 (bio) or asymmetric unit (asu).

        Returns:
            Path to the downloaded file.

        """
        pdb_id = pdb_id.upper()
        file_type = self._check_file_type(file_type)
        out_dir = (
            Path(self.out_dir / file_type / pdb_id[-3:-1])
            if self.subdir
            else Path(self.out_dir / file_type)
        )
        out_dir.mkdir(parents=True, exist_ok=True)
        gz_cif_file = out_dir / f"{pdb_id}.cif.gz"
        if gz_cif_file.exists():
            return gz_cif_file

        cif_url = (
            f"https://files.rcsb.org/download/{pdb_id}-assembly1.cif.gz"
            if file_type == "bio"
            else f"https://files.rcsb.org/download/{pdb_id}.cif.gz"
        )
        r = self.session.get(cif_url, timeout=self.timeout)
        r.raise_for_status()
        with open(gz_cif_file, "wb") as f:
            f.write(r.content)
        return gz_cif_file

    @staticmethod
    def _check_file_type(file_type: str) -> str:
        """Check if the file type is valid."""
        match file_type.lower():
            case "bio" | "biological" | "biological_assembly" | "biological-assembly":
                return "bio"
            case "asu" | "asymmetric" | "asymmetric_unit" | "asymmetric-unit":
                return "asu"
            case _:
                raise ValueError(f"Invalid file type: {file_type}")
