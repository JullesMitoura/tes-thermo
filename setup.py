from setuptools import setup
from setuptools import find_packages

with open("README.md", "r", encoding="utf-8") as arq:
    readme = arq.read()

setup(
    name='tes_thermo',
    version='0.1.7',
    license='MIT License',
    author='Julles Mitoura, Antonio Freitas and Adriano Mariano',
    author_email='mitoura96@outlook.com',
    description='TeS is a tool for simulating reaction processes. It uses the Gibbs energy minimization approach with the help of Pyomo and Ipopt as solvers.',
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
        'numpy',
        'scipy==1.16.0',
        'pyomo==6.9.2',
        'thermo==0.4.2',
        'faiss-cpu==1.11.0',
        'PyMuPDF==1.26.1',
        'langchain==0.3.13',
        'langchain_openai==0.2.14',
        'langchain_experimental==0.3.4'
    ],
)