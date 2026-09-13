#!/usr/bin/env python
# coding: utf-8

# In[48]:


var('n e_t u_bar v_bar f_4t f_5t g_4t g_5t Xi dt du dv gamma F_phi e_phi y K G m')
var('u_bar m G eta gamma e_t n Xi')
var('y y_0 h TE y_1PN y_15PN y_2PN y_25PN Xi chi_A chi_S delta omega')
var('e_t gamma y_inv e0 omega_dash A B flagTail e0 x0 flagHF SO chi w0 f0 j f_var f0_var e y0 kappaA kappaS horiz z')


# In[49]:


import re
filename = "henry/de_dt_final.txt"

with open(filename, "r") as f:
    content = f.read()


match = re.search(r"de_dt_final\s*=\s*(.*)", content, re.DOTALL)
if match:
    expr_str = match.group(1).strip()
    de_dt = ((sage_eval(expr_str, locals=globals())).subs({e:e_t})).expand()
    #*(1/(Xi^3*G))



filename = "henry/dx_dt_final.txt"

with open(filename, "r") as f:
    content = f.read()


match = re.search(r"dx_dt_final\s*=\s*(.*)", content, re.DOTALL)
if match:
    expr_str = match.group(1).strip()
    dx_dt = ((sage_eval(expr_str, locals=globals())).subs({e:e_t})).expand()
    #*(Xi^3/G)


filename = "henry/dx_dt_horiz.txt"

with open(filename, "r") as f:
    content = f.read()


match = re.search(r"dx_dt_horiz\s*=\s*(.*)", content, re.DOTALL)
if match:
    expr_str = match.group(1).strip()
    dx_dt_horiz = (((sage_eval(expr_str, locals=globals())).subs({e:e_t})).expand())*flagHF   


filename = "henry/de_dt_horiz.txt"

with open(filename, "r") as f:
    content = f.read()


match = re.search(r"de_dt_horiz\s*=\s*(.*)", content, re.DOTALL)
if match:
    expr_str = match.group(1).strip()
    de_dt_horiz = (((sage_eval(expr_str, locals=globals())).subs({e:e_t})).expand())*flagHF    


# In[50]:


filename = "henry/de_nonspinning_final.txt"

with open(filename, "r") as f:
    content = f.read()


match = re.search(r"de_nonspinning_final\s*=\s*(.*)", content, re.DOTALL)
if match:
    expr_str = match.group(1).strip()
    de_nonspinning_dt = (sage_eval(expr_str, locals=globals())).expand()

filename = "henry/dx_nonspinning_final.txt"

with open(filename, "r") as f:
    content = f.read()


match = re.search(r"dx_nonspinning_final\s*=\s*(.*)", content, re.DOTALL)
if match:
    expr_str = match.group(1).strip()
    dx_nonspinning_dt = (sage_eval(expr_str, locals=globals())).expand()


# In[51]:


import re
filename = "henry/xe.txt"

with open(filename, "r") as f:
    content = f.read()


match = re.search(r"xe_final\s*=\s*(.*)", content, re.DOTALL)
if match:
    expr_str = match.group(1).strip()
    xe_expr = ((sage_eval(expr_str, locals=globals()))*(e^(12/19)/e0^(12/19))).subs({e:e_t}).expand()


# In[52]:


def truncate_var(expr,var, order=2):
    """
    Truncate a Sage expression in Xi to the given order (default 3).
    Keeps terms up to Xi^order and discards higher powers.
    """
    truncated = sum(expr.coefficient(var,n) * var^n for n in range(order + 1))
    return truncated


# In[53]:


filename = "x_rederivation/x_2PN.txt"

with open(filename, "r") as f:
    content = f.read()

import re
match = re.search(r"x_2PN\s*=\s*(.*)", content, re.DOTALL)
if match:
    expr_str = match.group(1).strip()
    # Convert string back into a Sage expression
    x_2PN = sage_eval(expr_str, locals=globals())

x_N_1=x_2PN.coefficient(Xi,0)
x_1PN_final=x_2PN.coefficient(Xi,2)
x_15PN_final=x_2PN.coefficient(Xi,3)
x_2PN_final=x_2PN.coefficient(Xi,4)


# **Factor PN expressions to reduce ocmputational cost**

# **Factor the dx_dt expressions and include non spinning terms from Michael, Faye,...**

# In[54]:


dx_dt_N_term=(dx_dt.coefficient(Xi,0).subs({(-e_t^2+1)^(-7/2):y}).factor().subs({y:(-e_t^2+1)^(-7/2)}))


# In[55]:


dx_dt_1PN_term_base=(dx_dt.coefficient(Xi,2).coefficient(x^6).subs({(-e_t^2+1)^(-9/2):y}).factor().subs({y:(-e_t^2+1)^(-9/2)}))
dx_dt_1PN_term_num=(dx_dt_1PN_term_base.numerator().collect(e_t))*x^6
dx_dt_1PN_term_den=dx_dt_1PN_term_base.denominator()
dx_dt_1PN_term=dx_dt_1PN_term_num/dx_dt_1PN_term_den


# In[56]:


dx_dt_15PN_term_base=dx_dt.coefficient(Xi,3).coefficient(x^(13/2)).subs({e_t^10 - 5*e_t^8 + 10*e_t^6 - 10*e_t^4 + 5*e_t^2 - 1:-(-e_t^2+1)^5})
dx_dt_15PN_term_num=(dx_dt_15PN_term_base.numerator().collect(e_t))*x^(13/2)
dx_dt_15PN_term_den=dx_dt_15PN_term_base.denominator()
dx_dt_15PN_term=dx_dt_15PN_term_num/dx_dt_15PN_term_den


# In[57]:


dx_dt_2PN_term_base_0=(dx_nonspinning_dt.coefficient(Xi,4)+dx_dt.coefficient(Xi,4).coefficient(SO,1)*SO+dx_dt.coefficient(Xi,4).coefficient(SO,2)*SO^2)
dx_dt_2PN_term_base=(dx_dt_2PN_term_base_0.coefficient(x^7).subs({(-e_t^2+1):y}).subs({(e_t^2-1):-y}))
dx_dt_2PN_term_pt_1=(dx_dt_2PN_term_base.coefficient(y^(-5)).numerator().collect(e_t)/dx_dt_2PN_term_base.coefficient(y^(-5)).denominator())*(1/y^5)
dx_dt_2PN_term_pt_2=(dx_dt_2PN_term_base.coefficient(y^(-11/2)).numerator().collect(e_t)/dx_dt_2PN_term_base.coefficient(y^(-11/2)).denominator())*(1/y^(11/2))
dx_dt_2PN_term=((dx_dt_2PN_term_pt_1+dx_dt_2PN_term_pt_2).subs({y:-e_t^2+1}))*x^7


# In[58]:


dx_dt_25PN_term_base_0=(dx_dt.coefficient(Xi,5)+dx_dt_horiz.coefficient(Xi,5)).coefficient(x^(15/2)).subs({(-e_t^2+1)^(1/2):y^(1/2)}).subs({e_t^12 - 6*e_t^10 + 15*e_t^8 - 20*e_t^6 + 15*e_t^4 - 6*e_t^2 + 1:y^6})
dx_dt_25PN_term_pt_1=(dx_dt_25PN_term_base_0.coefficient(y^(-6)).numerator().collect(e_t)/dx_dt_25PN_term_base_0.coefficient(y^(-6)).denominator())*(1/y^6)
dx_dt_25PN_term_pt_2=(dx_dt_25PN_term_base_0.coefficient(y^(-11/2)).numerator().collect(e_t)/dx_dt_25PN_term_base_0.coefficient(y^(-11/2)).denominator())*(1/y^(11/2))


# In[59]:


dx_dt_25PN_term_base_1=(dx_nonspinning_dt.coefficient(Xi,5).coefficient(x^(15/2)).subs({-e_t^2+1:y}).subs({e_t^2-1:-y}))
dx_dt_25PN_term_pt_3=(dx_dt_25PN_term_base_1.subs({y:1/z}).coefficient(z,0).numerator().collect(e_t)/dx_dt_25PN_term_base_1.subs({y:1/z}).coefficient(z,0).denominator())
dx_dt_25PN_term_pt_4=(dx_dt_25PN_term_base_1.coefficient(y^(-1)).numerator().collect(e_t)/dx_dt_25PN_term_base_1.coefficient(y^(-1)).denominator())*(1/y)
dx_dt_25PN_term_pt_5=(dx_dt_25PN_term_base_1.coefficient(y^(-3/2)).numerator().collect(e_t)/dx_dt_25PN_term_base_1.coefficient(y^(-3/2)).denominator())*(1/y^(3/2))
dx_dt_25PN_term=((dx_dt_25PN_term_pt_1+dx_dt_25PN_term_pt_2+dx_dt_25PN_term_pt_3+dx_dt_25PN_term_pt_4+dx_dt_25PN_term_pt_5).subs({y:-e_t^2+1}))*x^(15/2)


# In[60]:


dx_dt_3PN_term_base_0=((dx_dt.coefficient(Xi,6)+dx_nonspinning_dt.coefficient(Xi,6)).expand().coefficient(x^8).subs({e_t^12 - 6*e_t^10 + 15*e_t^8 - 20*e_t^6 + 15*e_t^4 - 6*e_t^2 + 1:y^6}).subs({(-e_t^2+1):y}))
dx_dt_3PN_term_pt_1=(dx_dt_3PN_term_base_0.coefficient(y^(-6)).numerator().collect(e_t)/dx_dt_3PN_term_base_0.coefficient(y^(-6)).denominator())*(1/y^6)
dx_dt_3PN_term_pt_2=(dx_dt_3PN_term_base_0.coefficient(y^(-13/2)).numerator().collect(e_t)/dx_dt_3PN_term_base_0.coefficient(y^(-13/2)).denominator())*(1/y^(13/2))
dx_dt_3PN_term_pt_3=(dx_dt_3PN_term_base_0.subs({y:1/z}).coefficient(z,0).numerator().collect(e_t)/dx_dt_3PN_term_base_0.subs({y:1/z}).coefficient(z,0).denominator())
dx_dt_3PN_term=((dx_dt_3PN_term_pt_1+dx_dt_3PN_term_pt_2+dx_dt_3PN_term_pt_3).subs({y:-e_t^2+1}))*x^8


# In[61]:


dx_dt_final=dx_dt_N_term+dx_dt_1PN_term*Xi^2+dx_dt_15PN_term*Xi^3+dx_dt_2PN_term*Xi^4+dx_dt_25PN_term*Xi^5+dx_dt_3PN_term*Xi^6


# **Factor the de_dt expressions and include non spinning terms from Michael, Faye,...**

# In[62]:


de_dt_N_term=(de_dt.coefficient(Xi,0).subs({(-e_t^2+1)^(-5/2):y}).factor().subs({y:(-e_t^2+1)^(-5/2)}))


# In[63]:


de_dt_1PN_term_base=(de_dt.coefficient(Xi,2).coefficient(x^5).subs({(-e_t^2+1)^(-7/2):y}).factor().subs({y:(-e_t^2+1)^(-7/2)}))
de_dt_1PN_term_num=(de_dt_1PN_term_base.numerator().collect(e_t))*x^5
de_dt_1PN_term_den=de_dt_1PN_term_base.denominator()
de_dt_1PN_term=de_dt_1PN_term_num/de_dt_1PN_term_den


# In[64]:


de_dt_15PN_term_base=de_dt.coefficient(Xi,3).coefficient(x^(11/2)).subs({e_t^8 - 4*e_t^6 + 6*e_t^4 - 4*e_t^2 + 1:(-e_t^2+1)^4})
de_dt_15PN_term_num=(de_dt_15PN_term_base.numerator().collect(e_t))*x^(11/2)
de_dt_15PN_term_den=de_dt_15PN_term_base.denominator()
de_dt_15PN_term=de_dt_15PN_term_num/de_dt_15PN_term_den


# In[65]:


de_dt_2PN_term_base=de_dt.coefficient(Xi,4).coefficient(x^6).subs({e_t^8 - 4*e_t^6 + 6*e_t^4 - 4*e_t^2 + 1:y^4}).subs({(-e_t^2+1):y})
de_dt_2PN_term_pt_1=(de_dt_2PN_term_base.coefficient(y^(-9/2)).numerator().collect(e_t)/de_dt_2PN_term_base.coefficient(y^(-9/2)).denominator())*(1/y^(9/2))
de_dt_2PN_term_pt_2=(de_dt_2PN_term_base.coefficient(y^(-4)).numerator().collect(e_t)/de_dt_2PN_term_base.coefficient(y^(-4)).denominator())*(1/y^(4))
de_dt_2PN_term=((de_dt_2PN_term_pt_1+de_dt_2PN_term_pt_2).subs({y:-e_t^2+1}))*x^6


# In[66]:


de_dt_25PN_term_base_0=(de_dt.coefficient(Xi,5)+de_dt_horiz.coefficient(Xi,5)).coefficient(x^(13/2)).subs({e_t^10 - 5*e_t^8 + 10*e_t^6 - 10*e_t^4 + 5*e_t^2 - 1:-y^5}).subs({(-e_t^2+1):y})
de_dt_25PN_term_pt_1=(de_dt_25PN_term_base_0.coefficient(y^(-5)).numerator().collect(e_t)/de_dt_25PN_term_base_0.coefficient(y^(-5)).denominator())*(1/y^5)
de_dt_25PN_term_pt_2=(de_dt_25PN_term_base_0.coefficient(y^(-9/2)).numerator().collect(e_t)/de_dt_25PN_term_base_0.coefficient(y^(-9/2)).denominator())*(1/y^(9/2))


# In[67]:


de_dt_25PN_term_base_1=(de_nonspinning_dt.coefficient(Xi,5).coefficient(x^(13/2)).subs({-e_t^2+1:y}).subs({e_t^2-1:-y}))
de_dt_25PN_term_pt_3=(de_dt_25PN_term_base_1.coefficient(y^(-1/2)).numerator().collect(e_t)/de_dt_25PN_term_base_1.coefficient(y^(-1/2)).denominator())*(1/y^(1/2))
de_dt_25PN_term_pt_4=(de_dt_25PN_term_base_1.coefficient(y^(1/2)).numerator().collect(e_t)/de_dt_25PN_term_base_1.coefficient(y^(1/2)).denominator())*(y^(1/2))
de_dt_25PN_term_pt_5=(de_dt_25PN_term_base_1.coefficient(y,0).numerator().collect(e_t)/de_dt_25PN_term_base_1.coefficient(y,0).denominator())
de_dt_25PN_term=((de_dt_25PN_term_pt_1+de_dt_25PN_term_pt_2+de_dt_25PN_term_pt_3+de_dt_25PN_term_pt_4+de_dt_25PN_term_pt_5).subs({y:-e_t^2+1}))*x^(13/2)


# In[68]:


de_dt_3PN_term_base_0=((de_dt.coefficient(Xi,6)+de_nonspinning_dt.coefficient(Xi,6)).expand().coefficient(x^7).subs({e_t^10 - 5*e_t^8 + 10*e_t^6 - 10*e_t^4 + 5*e_t^2 - 1:-y^5}).subs({(-e_t^2+1):y}).subs({(e_t^2-1):-y}))
de_dt_3PN_term_base_1=(de_dt_3PN_term_base_0-(de_dt_3PN_term_base_0.coefficient((sqrt(y)+1)^(-1)))*((sqrt(y)+1)^(-1))).expand()
de_dt_3PN_term_base_2=de_dt_3PN_term_base_1+(de_dt_3PN_term_base_0.coefficient((sqrt(y)+1)^(-1)).subs({y^(-11/2):y_0}))*((sqrt(y)+1)^(-1))


# In[69]:


de_dt_3PN_term_pt_1=(de_dt_3PN_term_base_2.coefficient(y^(-5)).numerator().collect(e_t)/de_dt_3PN_term_base_2.coefficient(y^(-5)).denominator())*(1/y^5)
de_dt_3PN_term_pt_2=((de_dt_3PN_term_base_2.coefficient(y^(-11/2)).numerator().collect(e_t)/de_dt_3PN_term_base_2.coefficient(y^(-11/2)).denominator())+de_dt_3PN_term_base_2.coefficient(y_0))*(1/y^(11/2))
de_dt_3PN_term_pt_3=(de_dt_3PN_term_base_2.coefficient(y^(1/2)).numerator().collect(e_t)/de_dt_3PN_term_base_2.coefficient(y^(1/2)).denominator())*(y^(1/2))


# In[70]:


messy_term=(de_dt_3PN_term_base_2.coefficient(y,0)-(de_dt_3PN_term_base_2.coefficient(y_0))*(y_0)).expand()


# In[71]:


de_dt_3PN_term_pt_4=(messy_term.coefficient(y,0).numerator().collect(e_t)/messy_term.coefficient(y,0).denominator())
de_dt_3PN_term=((de_dt_3PN_term_pt_1+de_dt_3PN_term_pt_2+de_dt_3PN_term_pt_3+de_dt_3PN_term_pt_4).subs({y:-e_t^2+1}))*x^7


# In[72]:


de_dt_final=de_dt_N_term+de_dt_1PN_term*Xi^2+de_dt_15PN_term*Xi^3+de_dt_2PN_term*Xi^4+de_dt_25PN_term*Xi^5+de_dt_3PN_term*Xi^6


# **1 PN**

# In[73]:


dx_de_1PN_integrand=(((dx_dt_N_term+dx_dt_1PN_term*Xi^2)/(de_dt_N_term+de_dt_1PN_term*Xi^2))/x)


# In[74]:


dx_de_1PN_integrand_2=dx_de_1PN_integrand.taylor((Xi,0),(e_t,0),8).coefficient(Xi,2)


# **1.5 PN**

# In[77]:


dx_de_15PN_integrand=(((dx_dt_N_term+dx_dt_1PN_term*Xi^2+dx_dt_15PN_term*Xi^3)/(de_dt_N_term+de_dt_1PN_term*Xi^2+de_dt_15PN_term*Xi^3))/x)


# In[78]:


dx_de_15PN_integrand_2=dx_de_15PN_integrand.taylor((Xi,0),(e_t,0),9).coefficient(Xi,3)


# **2 PN**

# In[80]:


dx_de_2PN_integrand=(((dx_dt_N_term+dx_dt_1PN_term*Xi^2+dx_dt_15PN_term*Xi^3+dx_dt_2PN_term*Xi^4)/(de_dt_N_term+de_dt_1PN_term*Xi^2+de_dt_15PN_term*Xi^3+de_dt_2PN_term*Xi^4))/x)


# In[81]:


dx_de_2PN_integrand_2=dx_de_2PN_integrand.taylor((Xi,0),(e_t,0),10).coefficient(Xi,4)


# **Contributions**

# In[ ]:


dx_de_1_contribution=(dx_de_1PN_integrand_2.collect(x).subs({x:x_N_1+Xi^2*x_1PN_final+Xi^3*x_15PN_final+Xi^4*x_2PN_final}).taylor((Xi,0),(e_t,0),(e0,0),14))*Xi^2


# In[ ]:


filename = "x_rederivation/slow/dx_de_1_contribution.txt"


with open(filename, "w") as f:
    f.write("# Expansion of dx_de_1_contribution\n")
    f.write("# Plain ASCII output from Sage\n\n")
    f.write("dx_de_1_contribution =\n")
    f.write(str(dx_de_1_contribution))
    f.write("\n\n")


# In[ ]:


dx_de_15_contribution=(dx_de_15PN_integrand_2.collect(x).subs({x:x_N_1+Xi^2*x_1PN_final+Xi^3*x_15PN_final}).taylor((Xi,0),(e_t,0),(e0,0),12))*Xi^3


# In[ ]:


filename = "x_rederivation/slow/dx_de_15_contribution.txt"


with open(filename, "w") as f:
    f.write("# Expansion of dx_de_15_contribution\n")
    f.write("# Plain ASCII output from Sage\n\n")
    f.write("dx_de_15_contribution =\n")
    f.write(str(dx_de_15_contribution))
    f.write("\n\n")


# In[ ]:


dx_de_2_contribution=(dx_de_2PN_integrand_2.collect(x).subs({x:x_N_1+Xi^2*x_1PN_final}).taylor((Xi,0),(e_t,0),(e0,0),10))*Xi^4


# In[ ]:


filename = "x_rederivation/slow/dx_de_2_contribution.txt"


with open(filename, "w") as f:
    f.write("# Expansion of dx_de_2_contribution\n")
    f.write("# Plain ASCII output from Sage\n\n")
    f.write("dx_de_2_contribution =\n")
    f.write(str(dx_de_2_contribution))
    f.write("\n\n")


# In[ ]:




