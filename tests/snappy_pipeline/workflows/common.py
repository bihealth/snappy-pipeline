"""Shared method used in methods."""

import typing


def get_expected_log_files_dict(
    *, base_out: str, infix: typing.Optional[str] = None, extended: bool = False
):
    """
    :param base_out: Base path structure for log files. For example, if the expected path for
    the log is 'work/step.path/log/step.conda_info.txt', the argument should be
    'work/step.path/log/step'.
    :param infix: Optional infix string.
    :param extended: Whether to include env_yaml and wrapper in output

    :return: Returns dictionary with expected path for log files based on the provided input.
    """
    if infix:
        infix_dot = f"{infix}."
    else:
        infix_dot = ""
    expected = {
        "conda_info": f"{base_out}.{infix_dot}conda_info.txt",
        "conda_info_md5": f"{base_out}.{infix_dot}conda_info.txt.md5",
        "conda_list": f"{base_out}.{infix_dot}conda_list.txt",
        "conda_list_md5": f"{base_out}.{infix_dot}conda_list.txt.md5",
        "log": f"{base_out}.{infix_dot}log",
        "log_md5": f"{base_out}.{infix_dot}log.md5",
    }
    if extended:
        expected.update(
            {
                "env_yaml": f"{base_out}.{infix_dot}environment.yaml",
                "env_yaml_md5": f"{base_out}.{infix_dot}environment.yaml.md5",
                "wrapper": f"{base_out}.{infix_dot}wrapper.py",
                "wrapper_md5": f"{base_out}.{infix_dot}wrapper.py.md5",
            }
        )
    return expected


def get_expected_output_vcf_files_dict(base_out):
    """
    :param base_out: Base path structure for vcf files. For example, if the expected path for
    the log is 'work/step.path/log/step.vcf.gz', the argument should be
    'work/step.path/log/step'.
    :type base_out: str

    :return: Returns dictionary with expected path for vcf related files based on the
    provided input.
    """
    # Define expected
    expected = {
        "vcf": f"{base_out}.vcf.gz",
        "vcf_md5": f"{base_out}.vcf.gz.md5",
        "tbi": f"{base_out}.vcf.gz.tbi",
        "tbi_md5": f"{base_out}.vcf.gz.tbi.md5",
    }
    # Return
    return expected


def get_expected_output_bcf_files_dict(base_out):
    """
    :param base_out: Base path structure for log files. For example, if the expected path for
    the log is 'work/step.path/log/step.bcf', the argument should be
    'work/step.path/log/step'.
    :type base_out: str

    :return: Returns dictionary with expected path for bcf related files based on the
    provided input.
    """
    # Define expected
    expected = {
        "bcf": f"{base_out}.bcf",
        "bcf_md5": f"{base_out}.bcf.md5",
        "csi": f"{base_out}.bcf.csi",
        "csi_md5": f"{base_out}.bcf.csi.md5",
    }
    # Return
    return expected


def get_expected_gcnv_log_file(step_name):
    """
    :param step_name: Step name.
    :type step_name: str

    :return: Returns expected log file path for basic steps in gCNV.
    """
    expected_log = (
        "work/{mapper}.gcnv_"
        + step_name
        + ".{library_kit}/log/{mapper}.gcnv_"
        + step_name
        + ".{library_kit}.log"
    )
    return expected_log
