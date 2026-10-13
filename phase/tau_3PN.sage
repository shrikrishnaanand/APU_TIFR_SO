#!/usr/bin/env python
# coding: utf-8

# In[ ]:


var('n e_t u_bar v_bar f_4t f_5t g_4t g_5t Xi dt du dv gamma F_phi e_phi y K G m')
var('u_bar m G eta gamma e_t n Xi')
var('y y_0 h TE y_1PN y_15PN y_2PN y_25PN Xi chi_A chi_S delta omega')
var('e_t gamma y_inv e0 omega_dash A B flagTail e0 x0 flagHF SO chi w0 f0 j f_var f0_var e y0 kappaA kappaS horiz z')


# In[ ]:


import re
filename = "../x_e_files/x_3PN.txt"

with open(filename, "r") as f:
    content = f.read()


match = re.search(r"x_3PN\s*=\s*(.*)", content, re.DOTALL)
if match:
    expr_str = match.group(1).strip()
    xe_expr = (((sage_eval(expr_str, locals=globals()))).expand()).subs({SO^2:0},{SO^3:0})


# In[ ]:


import re
filename = "../henry/de_dt_final.txt"

with open(filename, "r") as f:
    content = f.read()


match = re.search(r"de_dt_final\s*=\s*(.*)", content, re.DOTALL)
if match:
    expr_str = match.group(1).strip()
    de_dt = ((sage_eval(expr_str, locals=globals())).subs({e:e_t})).expand()
    #*(1/(Xi^3*G))



filename = "../henry/dx_dt_final.txt"

with open(filename, "r") as f:
    content = f.read()


match = re.search(r"dx_dt_final\s*=\s*(.*)", content, re.DOTALL)
if match:
    expr_str = match.group(1).strip()
    dx_dt = ((sage_eval(expr_str, locals=globals())).subs({e:e_t})).expand()
    #*(Xi^3/G)


filename = "../henry/dx_dt_horiz.txt"

with open(filename, "r") as f:
    content = f.read()


match = re.search(r"dx_dt_horiz\s*=\s*(.*)", content, re.DOTALL)
if match:
    expr_str = match.group(1).strip()
    dx_dt_horiz = (((sage_eval(expr_str, locals=globals())).subs({e:e_t})).expand())*flagHF  


filename = "../henry/de_dt_horiz.txt"

with open(filename, "r") as f:
    content = f.read()


match = re.search(r"de_dt_horiz\s*=\s*(.*)", content, re.DOTALL)
if match:
    expr_str = match.group(1).strip()
    de_dt_horiz = (((sage_eval(expr_str, locals=globals())).subs({e:e_t})).expand())*flagHF   


# In[ ]:


filename = "../henry/de_nonspinning_final.txt"

with open(filename, "r") as f:
    content = f.read()


match = re.search(r"de_nonspinning_final\s*=\s*(.*)", content, re.DOTALL)
if match:
    expr_str = match.group(1).strip()
    de_nonspinning_dt = (sage_eval(expr_str, locals=globals())).expand()

filename = "../henry/dx_nonspinning_final.txt"

with open(filename, "r") as f:
    content = f.read()


match = re.search(r"dx_nonspinning_final\s*=\s*(.*)", content, re.DOTALL)
if match:
    expr_str = match.group(1).strip()
    dx_nonspinning_dt = (sage_eval(expr_str, locals=globals())).expand()


# In[ ]:


filename = "NNLO_e0/et_e0_NNLO_3PN.txt"

with open(filename, "r") as f:
    content = f.read()

# Extract the part after "dt_du_final ="
import re
match = re.search(r"et_e0_NNLO_3PN\s*=\s*(.*)", content, re.DOTALL)
if match:
    expr_str = match.group(1).strip()
    # Convert string back into a Sage expression
    et_e0_NNLO_3PN = (sage_eval(expr_str, locals=globals())).subs({SO^2:0},{SO^3:0})


# In[ ]:


dx_dt_horiz_term1=dx_dt_horiz.expand()
dx_dt_horiz_factor1=(dx_dt_horiz_term1.coefficient(chi_A^3)*chi_A^3+dx_dt_horiz_term1.coefficient(chi_S^3)*chi_S^3)


# In[ ]:


#This is still SO^3 term
dx_dt_horiz_term2=(dx_dt_horiz_term1-dx_dt_horiz_factor1).expand()
dx_dt_horiz_factor2=(dx_dt_horiz_term2.coefficient(chi_S^2)*chi_S^2+dx_dt_horiz_term2.coefficient(chi_A^2)*chi_A^2)


# In[ ]:


dx_dt_horiz_term3=(dx_dt_horiz_term2-dx_dt_horiz_factor2).expand()
dx_dt_horiz_factor3=(dx_dt_horiz_term3.coefficient(chi_S)*chi_S+dx_dt_horiz_term3.coefficient(chi_A)*chi_A)


# In[ ]:


dx_dt_horiz_labelled=(dx_dt_horiz_factor1*SO^3+dx_dt_horiz_factor2*SO^3+dx_dt_horiz_factor3*SO)


# In[ ]:


de_dt_horiz_term1=de_dt_horiz.expand()
de_dt_horiz_factor1=(de_dt_horiz_term1.coefficient(chi_A^3)*chi_A^3+de_dt_horiz_term1.coefficient(chi_S^3)*chi_S^3)


# In[ ]:


#This is still SO^3 term
de_dt_horiz_term2=(de_dt_horiz_term1-de_dt_horiz_factor1).expand()
de_dt_horiz_factor2=(de_dt_horiz_term2.coefficient(chi_S^2)*chi_S^2+de_dt_horiz_term2.coefficient(chi_A^2)*chi_A^2)


# In[ ]:


de_dt_horiz_term3=(de_dt_horiz_term2-de_dt_horiz_factor2).expand()
de_dt_horiz_factor3=(de_dt_horiz_term3.coefficient(chi_S)*chi_S+de_dt_horiz_term3.coefficient(chi_A)*chi_A)


# In[ ]:


de_dt_horiz_labelled=(de_dt_horiz_factor1*SO^3+de_dt_horiz_factor2*SO^3+de_dt_horiz_factor3*SO)


# **Factorize for faster runs**

# In[ ]:


dx_dt_N_term=(dx_dt.coefficient(Xi,0).subs({(-e_t^2+1)^(-7/2):y}).factor().subs({y:(-e_t^2+1)^(-7/2)}))


# In[ ]:


dx_dt_1PN_term_base=(dx_dt.coefficient(Xi,2).coefficient(x^6).subs({(-e_t^2+1)^(-9/2):y}).factor().subs({y:(-e_t^2+1)^(-9/2)}))
dx_dt_1PN_term_num=(dx_dt_1PN_term_base.numerator().collect(e_t))*x^6
dx_dt_1PN_term_den=dx_dt_1PN_term_base.denominator()
dx_dt_1PN_term=dx_dt_1PN_term_num/dx_dt_1PN_term_den


# In[ ]:


dx_dt_15PN_term_base=dx_dt.coefficient(Xi,3).coefficient(x^(13/2)).subs({e_t^10 - 5*e_t^8 + 10*e_t^6 - 10*e_t^4 + 5*e_t^2 - 1:-(-e_t^2+1)^5})
dx_dt_15PN_term_num=(dx_dt_15PN_term_base.numerator().collect(e_t))*x^(13/2)
dx_dt_15PN_term_den=dx_dt_15PN_term_base.denominator()
dx_dt_15PN_term=dx_dt_15PN_term_num/dx_dt_15PN_term_den


# In[ ]:


dx_dt_2PN_term_base_0=(dx_nonspinning_dt.coefficient(Xi,4)+dx_dt.coefficient(Xi,4).coefficient(SO,1)*SO+dx_dt.coefficient(Xi,4).coefficient(SO,2)*SO^2)
dx_dt_2PN_term_base=(dx_dt_2PN_term_base_0.coefficient(x^7).subs({(-e_t^2+1):y}).subs({(e_t^2-1):-y}))
dx_dt_2PN_term_pt_1=(dx_dt_2PN_term_base.coefficient(y^(-5)).numerator().collect(e_t)/dx_dt_2PN_term_base.coefficient(y^(-5)).denominator())*(1/y^5)
dx_dt_2PN_term_pt_2=(dx_dt_2PN_term_base.coefficient(y^(-11/2)).numerator().collect(e_t)/dx_dt_2PN_term_base.coefficient(y^(-11/2)).denominator())*(1/y^(11/2))
dx_dt_2PN_term=((dx_dt_2PN_term_pt_1+dx_dt_2PN_term_pt_2).subs({y:-e_t^2+1}))*x^7


# In[ ]:


dx_dt_25PN_term_base_0=(dx_dt.coefficient(Xi,5)+dx_dt_horiz_labelled.coefficient(Xi,5)).coefficient(x^(15/2)).subs({(-e_t^2+1)^(1/2):y^(1/2)}).subs({e_t^12 - 6*e_t^10 + 15*e_t^8 - 20*e_t^6 + 15*e_t^4 - 6*e_t^2 + 1:y^6})
dx_dt_25PN_term_pt_1=(dx_dt_25PN_term_base_0.coefficient(y^(-6)).numerator().collect(e_t)/dx_dt_25PN_term_base_0.coefficient(y^(-6)).denominator())*(1/y^6)
dx_dt_25PN_term_pt_2=(dx_dt_25PN_term_base_0.coefficient(y^(-11/2)).numerator().collect(e_t)/dx_dt_25PN_term_base_0.coefficient(y^(-11/2)).denominator())*(1/y^(11/2))


# In[ ]:


dx_dt_25PN_term_base_1=(dx_nonspinning_dt.coefficient(Xi,5).coefficient(x^(15/2)).subs({-e_t^2+1:y}).subs({e_t^2-1:-y}))
dx_dt_25PN_term_pt_3=(dx_dt_25PN_term_base_1.subs({y:1/z}).coefficient(z,0).numerator().collect(e_t)/dx_dt_25PN_term_base_1.subs({y:1/z}).coefficient(z,0).denominator())
dx_dt_25PN_term_pt_4=(dx_dt_25PN_term_base_1.coefficient(y^(-1)).numerator().collect(e_t)/dx_dt_25PN_term_base_1.coefficient(y^(-1)).denominator())*(1/y)
dx_dt_25PN_term_pt_5=(dx_dt_25PN_term_base_1.coefficient(y^(-3/2)).numerator().collect(e_t)/dx_dt_25PN_term_base_1.coefficient(y^(-3/2)).denominator())*(1/y^(3/2))
dx_dt_25PN_term=((dx_dt_25PN_term_pt_1+dx_dt_25PN_term_pt_2+dx_dt_25PN_term_pt_3+dx_dt_25PN_term_pt_4+dx_dt_25PN_term_pt_5).subs({y:-e_t^2+1}))*x^(15/2)


# In[ ]:


dx_dt_3PN_term_base_0=((dx_dt.coefficient(Xi,6)+dx_nonspinning_dt.coefficient(Xi,6)).expand().coefficient(x^8).subs({e_t^12 - 6*e_t^10 + 15*e_t^8 - 20*e_t^6 + 15*e_t^4 - 6*e_t^2 + 1:y^6}).subs({(-e_t^2+1):y}))
dx_dt_3PN_term_pt_1=(dx_dt_3PN_term_base_0.coefficient(y^(-6)).numerator().collect(e_t)/dx_dt_3PN_term_base_0.coefficient(y^(-6)).denominator())*(1/y^6)
dx_dt_3PN_term_pt_2=(dx_dt_3PN_term_base_0.coefficient(y^(-13/2)).numerator().collect(e_t)/dx_dt_3PN_term_base_0.coefficient(y^(-13/2)).denominator())*(1/y^(13/2))
dx_dt_3PN_term_pt_3=(dx_dt_3PN_term_base_0.subs({y:1/z}).coefficient(z,0).numerator().collect(e_t)/dx_dt_3PN_term_base_0.subs({y:1/z}).coefficient(z,0).denominator())
dx_dt_3PN_term=((dx_dt_3PN_term_pt_1+dx_dt_3PN_term_pt_2+dx_dt_3PN_term_pt_3).subs({y:-e_t^2+1}))*x^8


# In[ ]:


dx_dt_final=dx_dt_N_term+dx_dt_1PN_term*Xi^2+dx_dt_15PN_term*Xi^3+dx_dt_2PN_term*Xi^4+dx_dt_25PN_term*Xi^5+dx_dt_3PN_term*Xi^6


# **tau_3PN**

# In[ ]:


tau_3PN=(xe_expr/dx_dt_final.subs({SO^2:0},{SO^3:0})).subs({e_t:et_e0_NNLO_3PN}).taylor((Xi,0),(e0,0),13)


# In[ ]:


filename = "../phase/tau_files/tau_3PN.txt"


with open(filename, "w") as f:
    f.write("# Expansion of tau_3PN\n")
    f.write("# Plain ASCII output from Sage\n\n")
    f.write("x_3PN =\n")
    f.write(str(tau_3PN))
    f.write("\n\n")

