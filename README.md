## DQ fit library
### Python class
This python class is based on RooFit and allows to perform binned and unbinned fits. An example of fit can be found in /tutorial:

- Generate a toy sample
  ```ruby
  python tutorial.py configFit.json --gen_tutorial
  ```
- Run the fit on the toy sample
  ```ruby
  python tutorial.py configFit.json --do_fit
  ```
For validation of the code use the reference:
- [Analysis note Jpsi & Psi2S](https://alice-notes.web.cern.ch/system/files/notes/analysis/1216/2022-10-26-AN_Psi2S_v3.pdf)
