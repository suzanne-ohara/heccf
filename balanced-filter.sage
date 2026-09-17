r"""
Code to screen for cases where a heavenly elliptic curve could be *not* balanced, specifically in the
case where E is defined over a cubic number field.

REFERENCE: "Heavenly elliptic curves over cubic fields," 
           To appear.

           "Heavenly elliptic curves over quadratic fields,"
            https://arxiv.org/pdf/2410.18389

AUTHORS: Suzanne O'Hara (Wesleyan University)

Comments welcome: seohara 'typical email symbol' wesleyan 'typical email punctuation' edu
"""

# The possible values of e for an elliptic curve over a degree 3 number field.
e_val = [1, 2, 3, 4, 6, 8, 9, 12, 18]

# We hard-code the construction of possible Tate-Oort numbers (aka "j-vectors") for practical use:

poss_jvec = []
for e in e_val:
    j1 = 0
    while (j1 < e/2):
        j2 = e - j1
        poss_jvec.append( (j1, j2) )
        j1 += 1
        
# We now compute all the possible ranges of eigenvalues from the Weil bound for primes of norm q = p^f over Q.

def possible_frob_trace(p,f):
    q = p ** f
    extreme_frob_trace = floor(2 * q.sqrt())
    return list(range(-extreme_frob_trace, extreme_frob_trace +1))

# Next we define a recursive function that computes the values for trace of eth powers (tau_e) recursively.

def trace_frob_power( trace_frob, pow, prime_norm ):
    r'''Suppose theta is a Frobenius element for a prime pp of
        norm prime_norm, and that P(T) is the integer characteristic polynomial
        for theta. Let alpha, beta be the complex roots of P(T), so that

        trace_frob = alpha + beta.

        This function returns the integral value of the trace of theta^pow,
        which is necessarily

        alpha^pow + beta^pow.

        The approach is to work recursively with symmetric polynomials.'''

    # we use two recursions. First, if pow is even, we can work with pow/2.
    if pow % 2 == 0:
        # We use the identity (a^{2m} + b^{2m}) = (a^m + b^m)^2 - 2(ab)^m
        half_pow = pow/2
        return trace_frob_power( trace_frob, half_pow, prime_norm )**2 - 2 * prime_norm ** half_pow
    # otherwise, use a simpler recursion
    elif pow > 1:
        # We use the identity a^m + b^m = (a^{m-1} + b^{m-1})(a + b) - ab(a^{m-2} + b^{m-2})
        return trace_frob_power( trace_frob, pow - 1, prime_norm ) * trace_frob - prime_norm * trace_frob_power( trace_frob, pow - 2, prime_norm)
    elif pow == 1:
        return trace_frob
    else:
        raise Exception("pow must be a positive integer.")

# Then we define a function to compute the possible values of tau_e - q^{j1} - q^{j2}

def tau_comp( tau_e, q, j,k ):
    return (tau_e - q**j - q**k)

# Our first list of possible primes ell should be found as factors of tau_e - q^{j1} - q^{j2} computed above.
# We start with p = 3 to give a list of initial possibilites, and we will eliminate further ell values by comparing results from p = 2,5,7,11

poss_ell = {}

#The following loop computes tau_e for each possible frobenius trace of (3^f)
#Then uses this to compute and factor every possible value of tau_comp.

for jvec in poss_jvec:
    e = jvec[0] + jvec[1]
    poss_ell[jvec] = []
    for inertial_degree in [1,2,3]:
        prime_norm = 3 ** inertial_degree
        for tau in possible_frob_trace(3, inertial_degree):
            x = tau_comp(trace_frob_power(tau, e, prime_norm), prime_norm, jvec[0],jvec[1])
            for t in x.factor():
                if t[0] not in poss_ell[jvec]:
                    poss_ell[jvec].append(t[0])

#The primes remaining at this step should be the only ones where an unbalanced curve over a cubic field may exist.
#We can reduce the number of primes by comparing factors we get running these computations for other small primes.

#For book-keeping, we can re-order the lists of possible ell's to be from smallest to largest.
for jvec in poss_jvec:
    poss_ell[jvec].sort()


#In order for (j1, j2) to be a valid Tate-Oort pair for the prime ell, it must be that ell divides tau_e - p^{f*j1} - p^{f*j2} for any p != ell.
#So if tau_e - 5^{f*j1} - 5^{f*j2} != 0 (mod ell) for any f in {1,2,3} we mave remove the jvec pair from consideration
reduced_ell5= {}
for jvec in poss_jvec:
    reduced_ell5[jvec] = []
    e = jvec[0] + jvec[1]
    for inertial_degree in [1,2,3]:
        prime_norm = 5 ** inertial_degree
        for tau in possible_frob_trace(5, inertial_degree):
            x = tau_comp(trace_frob_power(tau, e, prime_norm), prime_norm, jvec[0], jvec[1])
            for y in poss_ell[jvec]:
                if x%y == 0 and y not in reduced_ell5[jvec]: reduced_ell5[jvec].append(y)


#reduce the possible list of jvec again against p = 7
reduced_ell7= {}
for jvec in poss_jvec:
    reduced_ell7[jvec] = []
    e = jvec[0] + jvec[1]
    for inertial_degree in [1,2,3]:
        prime_norm = 7 ** inertial_degree
        for tau in possible_frob_trace(7, inertial_degree):
            x = tau_comp(trace_frob_power(tau, e, prime_norm), prime_norm, jvec[0], jvec[1])
            for y in reduced_ell5[jvec]:
                if x%y == 0 and y not in reduced_ell7[jvec]: reduced_ell7[jvec].append(y)

#reduce the possible list of jvec again against p = 11
reduced_ell11= {}
for jvec in poss_jvec:
    reduced_ell11[jvec] = []
    e = jvec[0] + jvec[1]
    for inertial_degree in [1,2,3]:
        prime_norm = 11 ** inertial_degree
        for tau in possible_frob_trace(11, inertial_degree):
            x = tau_comp(trace_frob_power(tau, e, prime_norm), prime_norm, jvec[0], jvec[1])
            for y in reduced_ell7[jvec]:
                if x%y == 0 and y not in reduced_ell11[jvec]: reduced_ell11[jvec].append(y)

#reduce the possible list of jvec again against p = 2
reduced_ell2= {}
for jvec in poss_jvec:
    reduced_ell2[jvec] = []
    e = jvec[0] + jvec[1]
    for inertial_degree in [1,2,3]:
        prime_norm = 2 ** inertial_degree
        for tau in possible_frob_trace(2, inertial_degree):
            x = tau_comp(trace_frob_power(tau, e, prime_norm), prime_norm, jvec[0], jvec[1])
            for y in reduced_ell11[jvec]:
                if x%y == 0 and y not in reduced_ell2[jvec]: reduced_ell2[jvec].append(y)


#Now we remove the primes 2,3,5,7,11 from the list
#This is done because the above congruence on tau_e can only be done for p != ell

reduced_ell = reduced_ell2

for jvec in poss_jvec:
    reduced_ell[jvec].sort()
    for i in {2,3,5,7,11}:
        if i in reduced_ell[jvec]: reduced_ell[jvec].remove(i)

#Now for each remaining prime ell we find pairs (i1,i2) of powers of the l-adic cyclotomic character, chi, on the diagonal.
#These pairs should have the property that i1+i2 \equiv 1 (mod l-1) 
##and that the trace of frobenius element are in the appropriate range for the Weil bound

var('x')

my_cases = {}

for jvec in poss_jvec:
    e_ind = jvec[0] + jvec[1]
    s = reduced_ell[jvec]
    #for simplicity in calling things later, we re-name the list of possible primes ell that survive at the pair jvec

    my_cases[jvec] = {}

    for l in s:
        I_val = solve_mod((e_ind*x)==jvec[0],l-1)
        #We must have that e*i1 \equiv j1 (mod l-1)
        #This returns a tuple for each solution and we obtain i1 from the first entry

        poss_pairs = []
        #sets up a list for the possible (i1,i2) pairs of powers on the diagonal

        for i in I_val:
            poss_pairs.append((i[0], 1-i[0]))

        L = len(poss_pairs)
        #Later will loop over this length instead of elements of the lit
        ##Becaues we index through both trace3_f and poss_pairs at the same time.

        trace3_1 = []
        trace3_2 = []
        trace3_3 = []
        
        for I in poss_pairs: 
            trace3_1.append((3**I[0] + 3**I[1])%l)
            trace3_2.append((9**I[0] + 9**I[1])%l)
            trace3_3.append((27**I[0] + 27**I[1])%l)
        #compute the trace mod ell for all frob(3^f) options

        reduced_pairs = []
        #starts a counter for surviving (i1, i2) pairs

        
        for j in list(range(0,L)):
                if trace3_1[j] in possible_frob_trace(3,1) and poss_pairs[j] not in reduced_pairs: reduced_pairs.append(poss_pairs[j])
                elif -trace3_1[j]%l in possible_frob_trace(3,1) and poss_pairs[j] not in reduced_pairs: reduced_pairs.append(poss_pairs[j])
            #Checks if j or -j are within Weil bound range, if so, collect the pair it came in reduced_paris as a subset of poss_pairs
        
        for j in list(range(0,L)):
                if trace3_2[j] in possible_frob_trace(3,2) and poss_pairs[j] not in reduced_pairs: reduced_pairs.append(poss_pairs[j])
                elif -trace3_2[j]%l in possible_frob_trace(3,2) and poss_pairs[j] not in reduced_pairs: reduced_pairs.append(poss_pairs[j])

        for j in list(range(0,L)):
                if trace3_3[j] in possible_frob_trace(3,3) and poss_pairs[j] not in reduced_pairs: reduced_pairs.append(poss_pairs[j])
                elif -trace3_3[j]%l in possible_frob_trace(3,3) and poss_pairs[j] not in reduced_pairs: reduced_pairs.append(poss_pairs[j])
        
        if reduced_pairs != []:
            my_cases[jvec][l] = reduced_pairs
        #if there are no surviving (i1,i2) pairs, throw out the prime ell
        #otherwise, we keep ell and index the surviving (i1,i2) pairs


#Now that we have determined some (i1,i2) pairs that work, they must work for all primes p.
#We can remove pairs if the violate the Weil bound for other frobenius elements.

saved_cases = my_cases #I'm going to re-define this dictionary in the next loop, so we save the answers here just in case.

p = 5 

new_cases = {}


for jvec in poss_jvec:
    new_cases[jvec] = {}
    for l in my_cases[jvec].keys():
        new_cases[jvec][l] =[] #inintialize list of remaining (i1, i2) pairs

        for pair in my_cases[jvec][l]:
            tracep_1 = (p**(pair[0]) + p**(pair[1]))%l
            tracep_2 = ((p ** 2)**(pair[0]) + (p ** 2)**(pair[1]))%l
            tracep_3 = ((p ** 3)**(pair[0]) + (p ** 3)**(pair[1]))%l

            if tracep_1 in possible_frob_trace(p,1) and pair not in new_cases[jvec][l]: new_cases[jvec][l].append(pair)
            elif (-tracep_1)%l in possible_frob_trace(p,1) and pair not in new_cases[jvec][l]: new_cases[jvec][l].append(pair)
            elif tracep_2 in possible_frob_trace(p,2) and pair not in new_cases[jvec][l]: new_cases[jvec][l].append(pair)
            elif (-tracep_2)%l in possible_frob_trace(p,2) and pair not in new_cases[jvec][l]: new_cases[jvec][l].append(pair)
            elif tracep_3 in possible_frob_trace(p,3) and pair not in new_cases[jvec][l]: new_cases[jvec][l].append(pair)
            elif (-tracep_3)%l in possible_frob_trace(p,3) and pair not in new_cases[jvec][l]: new_cases[jvec][l].append(pair)

            #checked all possible frob traces against the correct Weil bound for each (p,f) pair.
            #(i1,i2) added to the list my_cases[jvec][l] only if it passes at least one of the above checks.

for p in [2,5,7,11,13,17]:
    for jvec in poss_jvec:
        for l in my_cases[jvec].keys():
            for pair in my_cases[jvec][l][:]:
                tracep_1 = (p**(pair[0]) + p**(pair[1]))%l
                tracep_2 = ((p ** 2)**(pair[0]) + (p ** 2)**(pair[1]))%l
                tracep_3 = ((p ** 3)**(pair[0]) + (p ** 3)**(pair[1]))%l

                if tracep_1 not in possible_frob_trace(p,1) and (-tracep_1)%l not in possible_frob_trace(p,1) and tracep_2 not in possible_frob_trace(p,2) and (-tracep_2)%l not in possible_frob_trace(p,2) and tracep_3 not in possible_frob_trace(p,3) and (-tracep_3)%l not in possible_frob_trace(p,3):
                    my_cases[jvec][l].remove(pair)


final_cases = {}
for jvec in poss_jvec:
    if my_cases[jvec] != {}:
        final_cases[jvec] = []
        for l in my_cases[jvec].keys():
            if my_cases[jvec][l] != [] and l not in final_cases[jvec]:
                final_cases[jvec].append(l)

#Note that final_cases drops a key jvec if there are no possible ell > 11 left.
#It also removes any keys ell that may have remained with an empty list of (i1,i2) pairs from the previous loop
#
#The final_cases dictionary reports Tate-Oort numbers paired with possible ell > 11 where it may be that
#the TO pair admidts an unbalanced heavenly curve at ell.
#
#Report the results to the user:
            

print("The remaining unbalanced cases with ell > 11 are as follows. \n")

print("Tate-Oort Numbers  Possible ell > 11")
print("-----------------  -----------------")   
for jvec in final_cases.keys():
    print( str(jvec).rjust(13), "      ", final_cases[jvec] )

#If one wishes to know which (i1,i2) pairs remain for some fixed Tate-Oort pair and prime
#Then for the non-empty lists can be found by using the following loop
#
#for jvec in final_cases.keys():
#    for l in final_cases[jvec]:
#        print(my_cases[jvec][l])

     

