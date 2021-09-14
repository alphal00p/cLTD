## Python implementation
The python script [cLTD.py](./cLTD.py) contains the function `integrate_energies` that generates the cLTD expression in a fully automatic fashion. 
Its input consists of a topology described in terms of a list of loop lines, where a loop line is the collection of all propagators with the same loop momenta flowing through it.
These loop lines are encoded in terms of the combination of loop momenta flowing though them (signature) and the number of propagators it contains (multiplicity). 


For the case of two-loop Feynman diagrams, any topology can be expressed in terms of three loop lines. 
One choice of the flows of the loop momenta can be $k$, $l$, and $k-l$, or written in terms of signatures for the basis $(k, l)$: `[[1, 0], [0, 1], [1, -1]]`.


### Example:
#### 2-loop Pentabox
The two-loop pentabox topology can be denoted in terms of the above signature together with the specification of the number of propagators per loop line: `[3, 1, 4]`.

The input to `integrate_energies` is then:
```python 
signatures = [[1, 0], [1, -1], [0, 1]]
n_props = [1, 3, 4]

res = integrate_energies(n_props, signatures,
                        verbose=False,
                        name="PentaBox",
                        output_type="mathematica")
```

which yields the cLTD expression. In this example the output consists of a mathematica file `cLTD_2l_pentabox.m`.
Other output formats are `yaml`, `pickle` (binary), and `FORM`.

#### 2-loop Sunrise
For the simpler topology of the 2-loop sunrise diagram we have as input:
```python
signatures = [[1, 0], [1, -1], [0, 1]]
n_props = [1, 1, 1]

res = integrate_energies(n_props, signatures,
                         verbose=False,
                         name='Sunrise',
                         output_type='FORM')
```
and we obtain as `FORM` output:
```python
L F = 
  -1*invd0*num(ncmd(+1*k1+1*E1-1*p1), ncmd(+1*E2-1*p2))
  -1*invd1*num(ncmd(+1*E0-1*p0), ncmd(+1*E0-1*p0+1*E1+1*p1))
  -1*num(ncmd(+1*E0-1*p0), ncmd(+1*E2-1*p2, +1*E0-1*p0+1*E1+1*p1))
  -1*num(ncmd(+1*E0-1*p0, +1*k1+1*E1-1*p1), ncmd(+1*E2-1*p2))
;
```
and
```python
d0 = +1*E0+1*p0+1*E1-1*p1+1*E2-1*p2;
d1 = +1*E0-1*p0+1*E1+1*p1+1*E2+1*p2;
```
where `invd<i>= 1/d<i>`. 

The `FORM` expression `F` is valid for arbitrary numerators in $k^0$ and $l^0$, the expression corresponding to a specific numerator has to be unfolded using [`pf_numerator.frm`](./pf_numerator.frm), where `num(ncmd(x1), ncmd(x2))` represents $\text{N}_\mathcal{N}([x_1],[x_2])$, as described in Appendix B of the paper.