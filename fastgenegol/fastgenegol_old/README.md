
## Python-C link:

Build genlife scheme with xcode using ../fastgenegol.xcodeproj/
Then
```
ln -s ../Build/Products/Release/libgenelife.dylib libgenelife.dylib
```
and
```
ln -s ../Build/Products/Release/libgenelifed.dylib libgenelifed.dylib
```

**Careful:**  make sure that Xcode is setup to compile inot the Release directory.  (Could be set to Debug directory, in which cased the links would have to change.


Then
```
python genelife.py
```

#### Variation by NP for gathering stats:
display ndisp, run for nrun without display, repeat:
|---ndisp---|------------------nrun---------------|
e.g. ndisp=100, nrun = 1000

python genelifeN.py

## Deterministic vs. nondeterministic code:

`genelife.ipynb` uses deterministic code.

`genelifed.ipynb` uses non-deterministic code.

`genelifed.ipynb` may be switched to non-deterministic by changing line 356 of `subgenelife4d.c` (to call either `update_nondet()` or `update_det()`), recompiling, and restarting the python kernel. A future commit can set this as a configuration parameter...


## Older, pre-python C link:

Activity graph:
```
cc -o actgenegol actgenegol.c

acitivity.py actgenegol
```

## Dynamics:

default: display 200; iterate 1000; repeat

```
cc -o dispgenegol dispgenegol.c

display.py dispgenegol
```

## Jupyter notebook

In the current directory ruyn `jupyter notebook`, then in the browser window that pops up look for `genelife.ipynb` and open it.  Execute cells down to the *Animation* section, and watch the animation.



