"""Shared method used in methods."""


def get_expected_log_files_dict(base_out):
    """
    :param base_out: Base path structure for log files. For example, if the expected path for
    the log is 'work/step.path/log/step.conda_info.txt', the argument should be
    'work/step.path/log/step'.
    :type base_out: str

    :return: Returns dictionary with expected path for log files based on the provided input.
    """
    # Define expected
    expected = {
        "conda_info": base_out + ".conda_info.txt",
        "conda_info_md5": base_out + ".conda_info.txt.md5",
        "conda_list": base_out + ".conda_list.txt",
        "conda_list_md5": base_out + ".conda_list.txt.md5",
        "log": base_out + ".log",
        "log_md5": base_out + ".log.md5",
    }
    # Return
    return expected
