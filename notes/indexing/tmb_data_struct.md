# Designing a data structure to be passed under TMB to calculate rate matrix in C++
This document is based on [the rate matrix model spec](https://github.com/mac-theobio/McMasterPandemic/blob/tmb/notes/indexing/doc-output/indexing.md) written by Steve Walker. Here, we are describing a data structure that can be transferred from the R side to the C++ side under TMB so that rate matrix can be calculated from the data structure on the C++ side. 

Here are some criteria that we are considering while designing such a data structure:

 - The data structure **must** contain sufficient information (such as state vector, parameters, and the "formula") to calculate the rate matrix.
 - The data structure **must** be losslessly transferrable from R to C++ under TMB.
 - The data structure should has as little redundancy as possible. The size of the data structure should be as small as possible.
 - Calculating a rate matrix from the data structure should be computationally efficient.
## Data transfer under TMB framework
First, lets review what data can be trasnferred under TMB. We can transfer data of
 - basic types such as integer, float, string, boolean
 - compound types such as vector, matrix (including sparse matrix)
 - aggregate type such as structure, which contains multiple variables of above types
 
Suppose you want to transfer 
   -    Vectors: v1 (as data) and p1 (as parameters)
   -    Sparse matrix: M 
    -   Scalar s1, 
    -   Boolean b1

On the R side, you will need 

    dd <- MakeADFun(data = list(
                            state = c(v1),
                            ratemat = M,
                            dt = s,
                            do_hazard = b1,
                            ),
                     parameters = list(params=c(p1)))

On the C++ side, you will need

    DATA_VECTOR(state);
    DATA_SPARSE_MATRIX(ratemat);
    DATA_INTEGER(dt);          
    DATA_INTEGER(do_hazard);            // boolean->integer
    PARAMETER_VECTOR(params);

**NOTE:** 
 1. We could encapsulate multiple data of various types into a structure and pass it to the C++ side and retrieve it by DATA_STRUCT(...). In this draft, we only consider transferring multiple data separately.
 2. Instead of passing data vector v1 and parameter vector p1 separately, we could concatenate them into one vector and pass the vector as data.

## Transferrable data structure
### On the R side
Let *s* be the state vector, *p* be the parameter vector, *M* be the rate matrix, ***sp*** (or `state_param` vector in Steve's jargon) be the vector concatenated of *s* and *p*. Also let *n* be the number of non-zero elements in *M*.

 - Create a pair of integer vectors ***from*** and ***to***, each of length *n*, containing the row and column indices of non-zero elements of *M* respectively. They are added to these two vectors in either row-major or column-major or other user-defined order. Whichever order is chosen, we will create the rest of vectors in that order.
 - Create an integer vector ***count*** of length *n*, where ***count**<sub>i</sub>*=number of operands (or *factors* in Steve's document) in the formula used to calculate the value of the *i-th* non-zero element of *M*
 - Create an integer vector ***spi*** by concatenating *n* vectors, each of which consists of a series of `state_param_index` of `ratemat_struct[[i]]`.
 - Create an integer vector ***modifier*** by concatenating *n* vectors, each of which consists of a series of numbers (in binary form): 0**00**, 0**01**, 0**10**, **1**xx (meaning identity, complement, inverse, addition operation respectively) corresponding to the variable referred by `state_param_index` of `ratemat_struct[[i]]`. By checking bits (bitwise AND), we can translate a segment of ***modifier*** vector back to formula (see section "*On the C++ side*" below).
### What to transfer
We then transfer ***sp*** as data (we can also transfer *s* and *p* separately and combine them on the C++ side if transferring ***sp*** breaks TMB) together with the vectors created in the above 4 steps. In other words, the data strcture to transfer is a collection of vectors:
 - ***sp***
 - ***from*** 
 - ***to***
 - ***count***
 - ***spi***
 - ***modifier***
### On the C++ side
After retrieving the transferred vectors, we can calculate rate matrix as follows (pseudo c++):

    int start = 0;
    int n = count.size();
    for (i=0; i<n; i++) {
        int row = from[i]; 
        int col = to[i];
        M[row,col] = 0;
        Type prod = 1;
        for (j=start; j<start+count[i]; j++) {
            Type x = sp[spi[j]];
            if (modifier[j] & 0b100) {
                M[row,col] += prod;
                prod = 1;
            }
            if (modifier[j] & 0b001)
                x = 1-x;
            else if (modifier[j] & 0b010)
                x = 1/x;
            prod *= x;
        }
        M[row,col] += prod;
        start += count[i];
    }

## Further considerations
### Formula generalization
In above sections, we assume the formula is in the form of sum of products of *factors*. For example

         (Ia) * (beta0) * (1/N) * (Ca) +
         (Ip) * (beta0) * (1/N) * (Cp) +
         (Im) * (beta0) * (1/N) * (Cm) * (1-iso_m) +
         (Is) * (beta0) * (1/N) * (Cs) * (1-iso_m))
If we want to reformat the above calculation into a more efficient formula:

         (Ia) * (beta0) * (1/N) * 
         (Ca + Cp + (Cm) * (1-iso_m) + (Cs) * (1-iso_m))
or we need a more flexible form, like

         (2.5*Ia - 17/N) + (1/Ca + 8) / (3 - iso_m)
The current proposed data structure cannot handle these cases. There are two possible solutions that I can think of:

 1. Creating a tree to represent a math expression with +, -, *, /
    operations. We will need to serialize the tree on the R side, pass
    it and transform it into a formula (in a similar way as
    deserializing it back to a tree) on the C++ side.
    
 2. Passing the formula as string. We can then evaluate the passed math expression using a math expression pasrer such as [Boost.Spirit](https://www.boost.org/doc/libs/1_76_0/libs/spirit/doc/html/index.html) or [exprtk](http://partow.net/programming/exprtk/index.html) or [TinyExpr](https://github.com/codeplea/tinyexpr)
## Summary
This document proposes a data structure and the associated rate matrix calculation method. The data strcture is filled on the R side, and transferred to and used in the C++ side. In the proposal, we actually implement a simple parser for a class of restricted math expressions.
    

