from setuptools import setup

if __name__ == "__main__":
    setup(package_data={'sensitivity_calculator': [
          'data/ACTPWV/ACT_annual_50.45.err', 'data/VariablePWV/*.out', 'data/*.txt', 'input.yaml']})
