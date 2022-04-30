# single-cell
Protocol for Single-Cell transcriptomic Analysis of Cell Fate Transitions during Cellular Reprogramming

## Setup

Make sure you have [anaconda](https://docs.anaconda.com/anaconda/install/index.html) installed on your machine.
Clone this repo, open a terminal and change your directory to the `single-cell` folder. Follow the below steps for setup.

```bash
conda create -y --name singlecell python=3.8
```

```bash
conda activate singlecell
```

```bash
pip install -r requirements.txt
```

## Running the analysis

Once you are setup, extract the files available in the `/data` folder before running the below command.

```bash
python scrun.py
```

> Figures from running the script show up in the `figures` folder.
