import pytest

def pytest_addoption(parser):
    parser.addoption("--numValidations", action="store", default=-1, help="the number of validation test datasets to run (default= -1 (all))")

@pytest.fixture(scope='session')
def numValidations(request):
    value = request.config.option.numValidations
    if value is None:
        value = -1
    return int(value)
