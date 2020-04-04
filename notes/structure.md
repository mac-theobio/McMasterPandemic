## notes on adding model structure

specifically, age and geography; also Erlang-y stuff

### general concerns

* adding stuff will slow the model down a lot. Maybe not a concern for forecasting with fixed params, but matters for estimation and parameter ensemble forecasting
* want to add structure in such a way that everything is still as transparent as possible, and that doesn't impede porting to a faster platform if necessary (see previous point)
    * allowing for parameters to vary appropriately across classes; make parameters into a list that can include vectors? (use unlist/relist)
	* keep capability to aggregate the big state space appropriately after running the sim (the current method can probably be extended, with some care - use structured names for subcompartments, e.g. Ia_a1_s2 for asymptomatic I in age class 1, spatial patch 2?)
	* appropriate use of Kronecker products to duplicate/structure the transition matrix (other possible stacking/unstacking tools? Probably _not_ Khatri-Rao, but vec-perm a la Caswell??)

### age

* age-structured data shouldn't be too hard
* WAIFW matrices from Prem et al 2013 if we want them
    * allows possibility of comparing control measures (school closure/work-from-home etc.), as in recent MRC paper
* age-specific mortality etc. from Riou?


### space

* need to implement spatial contact matrix (traffic flow?)
* age X space WAIFW can probably also work by Kronecker product (using "other" WAIFW matrix)
* county-level data also possible (county-by-age harder?)

### Erlang-y

* figure out convolutions?
* figure out how to auto-Erlangize a transition matrix (and then aggregate states); could allow for generalized Erlang/linear chain (Dushoff, Hurtado?)
