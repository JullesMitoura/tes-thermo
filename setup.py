from setuptools import setup
from setuptools import find_packages

with open("README.md", "r", encoding="utf-8") as arq:
    readme = arq.read()

setup(
    name='tes_thermo',
    version='0.1.9',
    license='MIT License',
    author='Julles Mitoura, Antonio Freitas and Adriano Mariano',
    author_email='mitoura96@outlook.com',
    description='TeS is a tool for simulating reaction processes. It uses the Gibbs energy minimization approach written as a nonlinear programming problem with Pyomo and is solved using IPOPT.',
    long_description=readme,
    long_description_content_type="text/markdown",
    keywords='gibbs, thermodynamics, virial, reactions, simulation, pyomo',
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "tes_thermo.solver": ["*.*", "**/*.*"],
    },
    install_requires=[
        'pandas==2.3.1',
        'numpy==2.4.0',
        'scipy==1.16.0',
        'pyomo==6.9.2',
        'thermo==0.6.0',
        'chemicals==1.5.0',
        'faiss-cpu==1.13.2',
        'PyMuPDF==1.26.1',
        'python-dotenv==1.2.1',
        'openai==2.14.0',
    ],
)