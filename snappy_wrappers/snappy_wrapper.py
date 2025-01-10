"""Abstract wrapper classes as utilities for snappy specific wrappers."""

import enum
import os
import shutil
import stat
import tempfile
import textwrap
from abc import abstractmethod, ABCMeta

from snakemake.shell import shell

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"


class SnappyWrapper(metaclass=ABCMeta):
    header = r"""
        # Pipe everything to log file
        if [[ -n "{snakemake.log.log}" ]]; then
            if [[ "$(set +e; tty; set -e)" != "" ]]; then
                rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
                exec &> >(tee -a "{snakemake.log.log}" >&2)
            else
                rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
                echo "No tty, logging disabled" >"{snakemake.log.log}"
            fi
        fi

        # Compute md5 except when filename ends with .md5
        compute_md5() {{
            fn=$1
            f=$(basename $fn)
            d=$(dirname $fn)
            pushd $d 1> /dev/null 2>&1
            md5sum $f > $f.md5
            popd 1> /dev/null 2>&1
        }}

        # Write out information about conda installation.
        conda list >{snakemake.log.conda_list}
        conda info >{snakemake.log.conda_info}
        compute_md5 {snakemake.log.conda_list}
        compute_md5 {snakemake.log.conda_info}

        # Create temp directory
        TMPDIR=$(mktemp -d)

        set -x

        # --------------------------------- Start command -----------------------------------------
    """

    footer = r"""
        # --------------------------------- End command -------------------------------------------

        set +x

        for fn in {snakemake.output}
        do
            if ! [[ $fn =~ \.md5$ || -p $fn ]]
            then
                compute_md5 $fn
            fi
        done
    """

    md5_log = r"""
        f=$(basename {log})
        d=$(dirname {log})
        pushd $d 1> /dev/null 2>&1
        md5sum $f > $f.md5
        popd 1> /dev/null 2>&1
    """

    output_links = r"""
        for path in {snakemake.output.output_links}; do
          dst=$path
          src=work/${{dst#output/}}
          if [[ ! -p $src ]]
          then
            ln -sr $src $dst
          fi
        done
    """

    def __init__(self, snakemake) -> None:
        self._snakemake = snakemake
        self._check_snakemake_attributes()

    def _check_snakemake_attributes(self):
        if not getattr(self._snakemake, "log", None):
            raise AttributeError("snakemake.log is not defined")
        if not getattr(self._snakemake.log, "log", None):
            raise AttributeError("snakemake.log.log is not defined")
        if not getattr(self._snakemake.log, "conda_list", None):
            raise AttributeError("snakemake.log.conda_list is not defined")
        if not getattr(self._snakemake.log, "conda_info", None):
            raise AttributeError("snakemake.log.conda_info is not defined")
        if not getattr(self._snakemake.log, "script", None):
            raise AttributeError("snakemake.log.script is not defined")

    @abstractmethod
    def run(self, cmd: str, do_md5: bool = True, with_output_links: bool = True) -> None:
        pass

    def _run(
        self, cmd: str, filename: str | None, do_md5: bool = True, with_output_links: bool = True
    ) -> None:
        """
        Creates a temp file for the script, executes it & computes the md5 sum of the log

        The shell script is first created as a temporary file, and then copied over to
        the log directory.
        This allows R scripts to be saved in the log directory, rather than the uninformative
        shell script starting R.

        : param cmd: The command string (after snakemake input/output/params expansion)
        : param filename: the path where to save the script
        """
        with tempfile.NamedTemporaryFile(mode="wt", delete_on_close=False) as f:
            tempfilename = f.name

            print(
                textwrap.dedent(
                    "\n".join(
                        (
                            SnappyWrapper.header.format(snakemake=self._snakemake),
                            cmd,
                            SnappyWrapper.footer.format(snakemake=self._snakemake)
                            if do_md5
                            else "",
                        )
                    )
                ),
                file=f,
            )

            f.flush()
            f.close()

            current_permissions = stat.S_IMODE(os.lstat(tempfilename).st_mode)
            os.chmod(tempfilename, current_permissions | stat.S_IXUSR)

            if filename is not None:
                shutil.copy(tempfilename, filename)

            shell(tempfilename)

        shell(SnappyWrapper.md5_log.format(log=str(self._snakemake.log.log)))

        if with_output_links and getattr(self._snakemake.output, "output_links", None) is not None:
            shell(SnappyWrapper.output_links.format(snakemake=self._snakemake))


class ShellWrapper(SnappyWrapper):
    def _run_bash(self, cmd: str, do_md5: bool = True, with_output_links: bool = True) -> None:
        self._run(
            cmd, self._snakemake.log.script, do_md5=do_md5, with_output_links=with_output_links
        )
        shell(SnappyWrapper.md5_log.format(log=self._snakemake.log.script))

    def run(self, cmd: str, do_md5: bool = True, with_output_links: bool = True) -> None:
        self._run_bash(cmd, do_md5=do_md5, with_output_links=with_output_links)


class RWrapper(SnappyWrapper):
    def _run_R(self, cmd: str, do_md5: bool = True, with_output_links: bool = True) -> None:
        with open(self._snakemake.log.script, "wt") as f:
            print(cmd, file=f)
        shell(SnappyWrapper.md5_log.format(log=self._snakemake.log.script))
        self._run(
            f"R --vanilla < {self._snakemake.log.script}",
            filename=None,
            do_md5=do_md5,
            with_output_links=with_output_links,
        )

    def run(self, cmd: str, do_md5: bool = True, with_output_links: bool = True) -> None:
        self._run_R(cmd, do_md5=do_md5, with_output_links=with_output_links)


# TODO: Put the BcftoolsWrapper (& BcftoolsCommand) in separate file


class BcftoolsCommand(enum.StrEnum):
    ANNOTATE = "annotate"
    CALL = "call"
    FILTER = "filter"
    MPILEUP = "mpileup"
    VIEW = "view"


class BcftoolsWrapper(ShellWrapper):
    additional_input_files = {
        BcftoolsCommand.ANNOTATE: (
            "annotations",
            "columns_file",
            "header_lines",
            "regions_file",
            "samples_file",
        ),
        BcftoolsCommand.CALL: (
            "ploidy_file",
            "regions_file",
            "samples_file",
            "targets_file",
            "group_samples",
        ),
        BcftoolsCommand.FILTER: ("mask_file", "regions_file", "targets_file"),
        BcftoolsCommand.MPILEUP: (
            "fasta_ref",
            "read_groups",
            "regions_file",
            "samples_file",
            "targets_file",
        ),
        BcftoolsCommand.VIEW: ("regions_file", "samples_file", "targets_file"),
    }

    bcftool_command = r"""
        bcftools {tool} \
            {extra_args} \
            {extra_files} \
            --output {snakemake.output.vcf} \
            {snakemake.input.vcf}
    """

    index_command = r"""
        if [[ ! -p {snakemake.output.vcf} ]]
        then
            tabix {snakemake.output.vcf}
        fi
    """

    def get_command(self, tool: BcftoolsCommand = BcftoolsCommand.VIEW):
        args = getattr(self._snakemake, "params", {}).get("args", {})
        extra_args = args.get("extra_args", [])
        extra_files = self._extra_input_files()
        cmd = BcftoolsWrapper.bcftool_command.format(
            tool=self._tool,
            snakemake=self._snakemake,
            extra_args=" \\\n    ".join(extra_args) if extra_args else "",
            extra_files=" \\\n    ".join(extra_files) if extra_files else "",
        )
        if args.get("index", False):
            cmd += "\n" + BcftoolsWrapper.index_command.format(snakemake=self._snakemake)
        return cmd

    def _extra_input_files(self) -> list[str]:
        extra_files = []
        for additional_input in BcftoolsWrapper.additional_input_files[self._tool]:
            filename = getattr(self._snakemake.input, additional_input, None)
            if filename is not None:
                extra_files.append(f"--{additional_input.replace('_', '-')} {filename}")
        return extra_files
