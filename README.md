# k2gap
A python module which contains some functions related to k2gap data. Most important being the selection function of targets selected by the K2GAP program (k2gap.sf). For a set of observables x, k2gap.sf returns the function y(x) which is True if it satisfies the selection criteria and False otherwise.



### Installation

    pip3 install  git+https://github.com/sanjibs/k2gap.git@main

### Example

    from k2gap import k2gap
    k2gap.sf(cno,ra,dec,jmag,hmag,kmag)
    help(k2gap)

OR

    import k2gap
    k2gap.sf(cno,ra,dec,jmag,hmag,kmag)
    help(k2gap.k2gap)

### Usage  

#### *Arguments*

* cno : (int) array_like, 
    - campaign number

* ra : (float) array_like, 
    - Right ascension [degree]

* dec : (float) array_like
    -  Declination  [degree]

* jmag : (float) array_like
    - 2MASS J band magnitude 

* hmag : (float) array_like
    - 2MASS H band magnitude 

* kmag : [float] array_like 
    - 2MASS Ks band magnitude 

* simulate_onsilicon : [boolean], optional keyword. 
    - If True it emulates the effect of K2fov.K2onSilicon.onSiliconCheck (https://github.com/KeplerGO/K2fov) by setting radius=1.4. The actual computation of onsilicion test is computationally expensive.  This is useful when working with mock data. See description >of keyword radius for more details.

* radius : [float] array_like, optional keyword
    - The distance of a star in [degree] from the center of its nearest 
    ccd module. The field of view of each campaign is divided up into 
    21 ccd modules with gaps between them. Circles defined by this radius and 
    centered on each module are non-overalping till a radius of 1.6 degrees.
    Use radius=1.75 to encapsulate the whole ccd module when working with 
    data for which K2fov.K2onSilicon.onSiliconCheck (https://github.com/KeplerGO/K2fov) 
    has been applied, e.g. observational data from K2. 
    Use radius=1.4, when working with 
    data for which K2fov.K2onSilicon.onSiliconCheck has not been applied, 
    e.g. a mock surveys. This has an area similar to that of a ccd module. 
    So using radius=1.4 is effectively similar to using radius=1.75 
    but with K2fov.K2onSilicon.onSiliconCheck.
    This is to avoid running K2fov.K2onSilicon.onSiliconCheck which is slow. 
    Note, for radius<=1.75 the four guiding ccds outside the science field of view 
    are not included in any of the circles.

#### *Returns*

* y : (boolean) array_like
    - 0/False if   

