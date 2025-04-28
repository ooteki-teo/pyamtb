from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="pyamtb",
    version="0.1.1",
    author="Wang Dinghui",
    author_email="wangdh@cumt.edu.cn",
    description="A Python package for tight-binding model calculations for altermagnets",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ooteki-teo/pyamtb",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.18.0",
        "matplotlib>=3.0.0",
        "tomlkit>=0.10.0",
        "pythtb>=1.7.2",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.9",
        ],
    },
    include_package_data=True,
    package_data={
        "pyamtb": ["*.toml"],
    },
    entry_points={
        "console_scripts": [
            "pyamtb=pyamtb.cli:main",
        ],
    },
) 