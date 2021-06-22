#!/usr/bin/env python
"""
Run HOMER's findMotifsGenome.pl program, or generate commands to do so.

Without -X/--execute, this script will simply produce (print) the command
needed to run HOMER's motif finding program for each sample. If that flag is
specified, though, these commands will be submitted to run.

Usage:
    /path/to/this/script [OPTIONS] --genome <assembly-name> <INPUT1> ...

Other useful options:
    1. -X/--execute: actually run commands, rather than just generating.
    2. --bsub: Path to LSF/bsub configuration file (resource specification
    for each command/job).
    3. --homer-path: Path to a specific HOMER installation; if unspecified,
    it's assumed that the appropriate program is on PATH
    (i.e., findMotifsGenome.pl will run).
    4. --output-prefix: Prefix for each HOMER directory created.

Outputs:
    Directly, the commands are printed. If run, for each input file, HOMER
    will create a folder with a corresponding name, with all of the
    motif finding output for that particular sample/input file.

"""

import argparse
from functools import partial
import os
import subprocess
import sys


__author__ = "Vince Reuter"
__email__ = "reuterv1@mskcc.org"



def _parse_cmdl(cmdl):
    """ Define and parse command-line interface. """

    parser = argparse.ArgumentParser(
        description="Run HOMER's findMotifsGenome.pl program",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required
    parser.add_argument(
        "input_files", nargs='+',
        help="Path(s) to input file(s), each of which will be used as target "
             "sequences for motif finding. These should be in BED or HOMER "
             "peak file format.")
    parser.add_argument(
        "-G", "--genome", required=True,
        help="Genomic assembly to used for the motif finding")

    # Optional
    parser.add_argument(
        "--bsub",
        help="Path to bsub configuration file; only use this if you want "
             "each command to be seen or submitted as a job on a cluster "
             "that uses LSF as the job scheduler.")
    parser.add_argument(
        "--homer-path",
        help="Explicit specification of path to HOMER installation; since "
             "the HOMER installation recommends making this part of PATH, "
             "this isn't a required option, as the HOMER programs are assumed "
             "to be available on PATH")
    parser.add_argument(
        "--logfile-name", default="motifs.log",
        help="Name used for logfile placed in each HOMER output folder; "
             "only relevant for cluster submissions")
    parser.add_argument(
        "-N", "--n-motifs", type=int, default=25,
        help="Number of motifs to find; HOMER docs recommend, if anything, "
             "reducing this default--especially for long motifs--to reduce "
             "runtime")
    parser.add_argument(
        "--output-prefix", default="Motifs_",
        help="Prefix for each output folder")
    parser.add_argument(
        "--region-size", default=200,
        help="Size of regions to use in the motif finding")
    parser.add_argument(
        "-X", "--execute", action="store_true",
        help="Actually execute commands, don't just print.")

    return parser.parse_args(cmdl)



def build_command(infile, outdir, genome, region_size, n_motifs, prog):
    """
    Build a command for the findMotifsGenome program.

    Parameters
    ----------
    infile : str
        Path to peak calls file.
    outdir : str
        Path to folder for output.
    genome : str
        Name of assembly to use for motif finding (should accord with
        what was used to call the peaks in the input file).
    region_size : int or str
        Number of base pairs in regions to search for motifs, or 'given' to
        simply use the size implied by the coordinates of the peak call.
    n_motifs : int
        Number of motifs to find.
    prog : str
        Path to the findMotifsGenome.pl program.

    Returns
    -------
    str
        Command with which to run HOMER's findMotifsGenome.pl program.

    """
    return "{p} {i} {g} {o} -size {s} -S {n}".format(
        p=prog, i=infile, o=outdir, g=genome, s=region_size, n=n_motifs)



def parse_bscub_cfg(bsub_cfg_file):
    """
    Parse a bsub configuration file, assigning values to submission params.

    Parameters
    ----------
    bsub_cfg_file : str
        Path to bsub configuration file

    Returns
    -------
    Mapping
        Value for each of a collection of LSF job submission parameters.

    """
    _, filetype = os.path.splitext(bsub_cfg_file)
    with open(bsub_cfg_file, 'r') as cfg:
        if filetype == ".json":
            import json
            conf_data = json.load(cfg)
        elif filetype == ".yaml":
            import yaml
            conf_data = yaml.load(cfg)
        elif filetype in [".csv", ".txt", ".tsv", ".tab"]:
            sep = "," if filetype == ".csv" else "\t"
            conf_data = {}
            for l in cfg:
                try:
                    k, v = l.strip().split(sep)
                except ValueError:
                    print("Skipping line: {}".format(l))
                else:
                    conf_data[k] = v
        else:
            raise ValueError("Unsupported filetype for bsub config: {}".
                             format(bsub_cfg_file))
        return conf_data



def main(cmdl):
    """ Run the script. """
    opts = _parse_cmdl(cmdl)
    try:
        region_size = int(opts.region_size)
    except (TypeError, ValueError):
        if opts.region_size != "given":
            raise ValueError(
                "Invalid size option to HOMER (should be an integer or "
                "'given'): {}".format(opts.region_size))
        region_size = opts.region_size
    prog = "findMotifsGenome.pl"
    if opts.homer_path:
        prog = os.path.join(opts.homer_path, prog)
    get_cmd = partial(
        build_command, genome=opts.genome, region_size=region_size,
        n_motifs=opts.n_motifs, prog=prog)
    if opts.bsub:
        print("Parsing LSF submission configurtation: {}".format(opts.bsub))
        bsub_conf = parse_bscub_cfg(opts.bsub) if opts.bsub else {}
        bsub_cmd_base = "bsub -n{slots} -W {time} -R \\\"rusage[mem={memory}]\\\"".\
                format(**bsub_conf)
    else:
        bsub_cmd_base = ""
    for infile in opts.input_files:
        basepath, _ = os.path.splitext(infile)
        input_folder, basename = os.path.split(basepath)
        outname = opts.output_prefix + basename
        outpath = os.path.join(input_folder, outname) \
            if input_folder else outname
        if not os.path.exists(outpath):
            os.makedirs(outpath)
        cmd = get_cmd(infile, outdir=outpath)
        if bsub_cmd_base:
            bsub_outfile = os.path.join(outpath, opts.logfile_name)
            cmd = "{} -o {} {}".format(bsub_cmd_base, bsub_outfile, cmd)
        if opts.execute:
            print("Submitting:\n{}".format(cmd))
            subprocess.check_call([part.replace('\\"', '\"') for part in cmd.split()])
        else:
            print("Dry run, no submission:\n{}\n".format(cmd))



if __name__ == "__main__":
    main(sys.argv[1:])
