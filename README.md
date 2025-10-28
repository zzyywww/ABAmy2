# ABAmy2
AB-Amy 2.0 offers a better and more robust tool for amyloidogenic risk prediction of therapeutic antibody light chain


you can clone this repository locally:
```
git clone https://github.com/zzyywww/ABAmy2.git 
cd ABAmy2
pip install -r requirements.txt --ignore-installed
```
Usage:

```
python ABAmy2.py [inputfile] [outputfile]
```

inputfile: **Fv region** of antibody light chain with **.fasta** format

outputfile: the predict result

Example:
```
python ABAmy2.py ./data/test.fasta ./result.txt
```

Paper describing this work has been published in Methods.

**Cited as**: Zhou Y, Liu W, Luo C, Huang Z, Samarappuli Mudiyanselage Savini G, Zhao L, Wang R, Huang J. Ab-Amy 2.0: Predicting light chain amyloidogenic risk of therapeutic antibodies based on antibody language model. Methods. 2025, 233:11-18. doi: 10.1016/j.ymeth.2024.11.005.
