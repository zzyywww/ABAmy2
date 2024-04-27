# ABAmy2
AB-Amy 2.0 offers a better and more robust tool for amyloidogenic risk prediction of therapeutic antibody light chain


you can clone this repository locally:
```
git clone https://github.com/zzyywww/ABAmy2.git 
cd ABAmy2
install -r requirements.txt --ignore-installed
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
