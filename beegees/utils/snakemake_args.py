"""Build the snakemake subprocess command for the BeeGees pipeline."""
from pathlib import Path

from beegees.utils.configs import get_snakefile, get_bundled_profile


def build_snakemake_cmd(
    configfile: Path,
    cores: int | None,
    profile: str | None,
    dryrun: bool,
    unlock: bool,
    extra_args: list[str],
) -> list[str]:
    cmd = ["snakemake", "--snakefile", str(get_snakefile())]

    cmd += ["--configfile", str(configfile)]

    if profile:
        bundled = get_bundled_profile(profile)
        cmd += ["--profile", str(bundled) if bundled.is_dir() else profile]
    else:
        cmd += ["--profile", str(get_bundled_profile("local"))]

    if cores:
        cmd += ["--cores", str(cores)]
    if dryrun:
        cmd += ["--dryrun"]
    if unlock:
        cmd += ["--unlock"]
    else:
        cmd += ["--rerun-incomplete"]

    return cmd + extra_args
