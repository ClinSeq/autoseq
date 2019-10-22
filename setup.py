from setuptools import setup, find_packages
from pip._internal.req import parse_requirements

# parse_requirements() returns generator of pip.req.InstallRequirement objects
install_reqs = parse_requirements("requirements.txt", session=False)

# reqs is a list of requirement
reqs = [str(ir.req) for ir in install_reqs if ir.req is not None]

setup(name='autoseq',
      version='1.0.0',
      packages=find_packages(exclude=('tests*', 'docs', 'examples')),
      install_requires=reqs,
      entry_points={
          'console_scripts': [
              'autoseq = autoseq.cli.cli:cli',
              'report2json = autoseq.report2json:main',
              'generate-ref = autoseq.generate_ref:main',
              'jobs2gantt = autoseq.cli.jobs2gantt:cli'
          ]
      }
      )
