import pytest

from snappy_wrappers.wrapper_parallel import (
    days,
    gib,
    gib_to_string,
    hours,
    kib,
    kib_to_string,
    mib,
    mib_to_string,
    minutes,
)

# Test isolated methods ----------------------------------------------------------------------------


def test_kib():
    """Tests wrapper_parallel.kib() call."""
    # Define expected dictionary
    expected_dict = {1: 1024, 2: 2048, 3: 3072, 1.5: 1024}
    # Get actual and assert
    for input_, expected in expected_dict.items():
        actual = kib(input_)
        assert actual == expected


def test_kib_exception():
    """Tests wrapper_parallel.kib() exception."""
    # Test raise ValueError
    with pytest.raises(ValueError):
        kib("One")
    # Test raise TypeError
    bad_input_list = [{1}, [1, 2]]
    for input_ in bad_input_list:
        with pytest.raises(TypeError):
            kib(input_)


def test_kib_to_string():
    """Tests wrapper_parallel.kib_to_string() call."""
    # Define expected dictionary
    expected_dict = {1: "1k", 2: "2k", 3: "3k", 1.5: "1k"}
    # Get actual and assert
    for input_, expected in expected_dict.items():
        actual = kib_to_string(input_)
        assert actual == expected


def test_mib():
    """Tests wrapper_parallel.mib() call."""
    # Define expected dictionary
    expected_dict = {1: 1048576, 2: 2097152, 3: 3145728, 1.5: 1048576}
    # Get actual and assert
    for input_, expected in expected_dict.items():
        actual = mib(input_)
        assert actual == expected


def test_mib_exception():
    """Tests wrapper_parallel.mib() exception."""
    # Test raise ValueError
    with pytest.raises(ValueError):
        mib("One")
    # Test raise TypeError
    bad_input_list = [{1}, [1, 2]]
    for input_ in bad_input_list:
        with pytest.raises(TypeError):
            mib(input_)


def test_mib_to_string():
    """Tests wrapper_parallel.mib_to_string() call."""
    # Define expected dictionary
    expected_dict = {1: "1M", 2: "2M", 3: "3M", 1.5: "1M"}
    # Get actual and assert
    for input_, expected in expected_dict.items():
        actual = mib_to_string(input_)
        assert actual == expected


def test_gib():
    """Tests wrapper_parallel.gib() call."""
    # Define expected dictionary
    expected_dict = {1: 1073741824, 2: 2147483648, 3: 3221225472, 1.5: 1073741824}
    # Get actual and assert
    for input_, expected in expected_dict.items():
        actual = gib(input_)
        assert actual == expected


def test_gib_exception():
    """Tests wrapper_parallel.gib() exception."""
    # Test raise ValueError
    with pytest.raises(ValueError):
        gib("One")
    # Test raise TypeError
    bad_input_list = [{1}, [1, 2]]
    for input_ in bad_input_list:
        with pytest.raises(TypeError):
            gib(input_)


def test_gib_to_string():
    """Tests wrapper_parallel.gib_to_string() call."""
    # Define expected dictionary
    expected_dict = {1: "1G", 2: "2G", 3: "3G", 1.5: "1G"}
    # Get actual and assert
    for input_, expected in expected_dict.items():
        actual = gib_to_string(input_)
        assert actual == expected


def test_minutes():
    """Tests wrapper_parallel.minutes()"""
    # Define expected values
    expected_dict = {1: "0:01:00", 0.5: "0:00:30"}
    # Get actual values and assert
    for input_, expected in expected_dict.items():
        actual = minutes(input_)
        assert actual == expected


def test_hours():
    """Tests wrapper_parallel.hours()"""
    # Define expected values
    expected_dict = {
        1: "1:00:00",
        0.5: "0:30:00",
        24: "1-00:00:00",
        25: "1-01:00:00",
        47: "1-23:00:00",
        48: "2-00:00:00",
    }
    # Get actual values and assert
    for input_, expected in expected_dict.items():
        actual = hours(input_)
        assert actual == expected


def test_days():
    """Tests wrapper_parallel.days()"""
    # Define expected values
    expected_dict = {1: "1 day, 0:00:00", 0.5: "12:00:00"}
    # Get actual values and assert
    for input_, expected in expected_dict.items():
        actual = days(input_)
        assert actual == expected
