#!/usr/bin/env python
# coding: utf-8

# In[1]:


var('n e_t u_bar v_bar f_4t f_5t g_4t g_5t Xi dt du dv gamma F_phi e_phi y K G m')
var('u_bar m G eta gamma e_t n Xi')
var('y y_0 h TE y_1PN y_15PN y_2PN y_25PN Xi chi_A chi_S delta omega')
var('e_t gamma y_inv e0 omega_dash A B flagTail e0 x0 flagHF SO chi w0 f0 j f_var f0_var e y0 kappaA kappaS horiz z')


# In[2]:


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


# In[3]:


dx_dt_horiz_term1=dx_dt_horiz.expand()
dx_dt_horiz_factor1=(dx_dt_horiz_term1.coefficient(chi_A^3)*chi_A^3+dx_dt_horiz_term1.coefficient(chi_S^3)*chi_S^3)


# In[4]:


#This is still SO^3 term
dx_dt_horiz_term2=(dx_dt_horiz_term1-dx_dt_horiz_factor1).expand()
dx_dt_horiz_factor2=(dx_dt_horiz_term2.coefficient(chi_S^2)*chi_S^2+dx_dt_horiz_term2.coefficient(chi_A^2)*chi_A^2)


# In[5]:


dx_dt_horiz_term3=(dx_dt_horiz_term2-dx_dt_horiz_factor2).expand()
dx_dt_horiz_factor3=(dx_dt_horiz_term3.coefficient(chi_S)*chi_S+dx_dt_horiz_term3.coefficient(chi_A)*chi_A)


# In[6]:


dx_dt_horiz_labelled=(dx_dt_horiz_factor1*SO^3+dx_dt_horiz_factor2*SO^3+dx_dt_horiz_factor3*SO)


# In[8]:


de_dt_horiz_term1=de_dt_horiz.expand()
de_dt_horiz_factor1=(de_dt_horiz_term1.coefficient(chi_A^3)*chi_A^3+de_dt_horiz_term1.coefficient(chi_S^3)*chi_S^3)


# In[9]:


#This is still SO^3 term
de_dt_horiz_term2=(de_dt_horiz_term1-de_dt_horiz_factor1).expand()
de_dt_horiz_factor2=(de_dt_horiz_term2.coefficient(chi_S^2)*chi_S^2+de_dt_horiz_term2.coefficient(chi_A^2)*chi_A^2)


# In[10]:


de_dt_horiz_term3=(de_dt_horiz_term2-de_dt_horiz_factor2).expand()
de_dt_horiz_factor3=(de_dt_horiz_term3.coefficient(chi_S)*chi_S+de_dt_horiz_term3.coefficient(chi_A)*chi_A)


# In[11]:


de_dt_horiz_labelled=(de_dt_horiz_factor1*SO^3+de_dt_horiz_factor2*SO^3+de_dt_horiz_factor3*SO)


# In[13]:


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


# In[14]:


import re
filename = "../henry/xe.txt"

with open(filename, "r") as f:
    content = f.read()


match = re.search(r"xe_final\s*=\s*(.*)", content, re.DOTALL)
if match:
    expr_str = match.group(1).strip()
    xe_expr = ((sage_eval(expr_str, locals=globals()))*(e^(12/19)/e0^(12/19))).subs({e:e_t}).expand()


# In[15]:


def truncate_var(expr,var, order=2):
    """
    Truncate a Sage expression in Xi to the given order (default 3).
    Keeps terms up to Xi^order and discards higher powers.
    """
    truncated = sum(expr.coefficient(var,n) * var^n for n in range(order + 1))
    return truncated


# In[16]:


#Labelling this x_2PN caused huge error. Be careful. Renamed to x_2PN_expr
filename = "../x_e_files/x_2PN.txt"

with open(filename, "r") as f:
    content = f.read()


match = re.search(r"x_2PN\s*=\s*(.*)", content, re.DOTALL)
if match:
    expr_str = match.group(1).strip()
    # Convert string back into a Sage expression
    x_2PN_expr = (sage_eval(expr_str, locals=globals())).expand()

x_N_1=x_2PN_expr.coefficient(Xi,0)
x_1PN_final=x_2PN_expr.coefficient(Xi,2)
x_15PN_final=x_2PN_expr.coefficient(Xi,3)
x_2PN_final=x_2PN_expr.coefficient(Xi,4)


# In[17]:


filename = "../x_e_files/slow/dx_de_1_contribution.txt"

with open(filename, "r") as f:
    content = f.read()

match = re.search(r"dx_de_1_contribution\s*=\s*(.*)", content, re.DOTALL)
if match:
    expr_str = match.group(1).strip()
    dx_de_1_contribution = ((sage_eval(expr_str, locals=globals())).expand())


# In[18]:


filename = "../x_e_files/slow/dx_de_15_contribution.txt"

with open(filename, "r") as f:
    content = f.read()
match = re.search(r"dx_de_15_contribution\s*=\s*(.*)", content, re.DOTALL)
if match:
    expr_str = match.group(1).strip()
    dx_de_15_contribution = ((sage_eval(expr_str, locals=globals())).expand())


# In[19]:


filename = "../x_e_files/slow/dx_de_2_contribution.txt"

with open(filename, "r") as f:
    content = f.read()

match = re.search(r"dx_de_2_contribution\s*=\s*(.*)", content, re.DOTALL)
if match:
    expr_str = match.group(1).strip()
    dx_de_2_contribution = ((sage_eval(expr_str, locals=globals())).expand())


# In[20]:


e_factor=(e_t^(12/19)/e0^(12/19))


# **Factorize for faster runs**

# In[21]:


dx_dt_N_term=(dx_dt.coefficient(Xi,0).subs({(-e_t^2+1)^(-7/2):y}).factor().subs({y:(-e_t^2+1)^(-7/2)}))


# In[22]:


dx_dt_1PN_term_base=(dx_dt.coefficient(Xi,2).coefficient(x^6).subs({(-e_t^2+1)^(-9/2):y}).factor().subs({y:(-e_t^2+1)^(-9/2)}))
dx_dt_1PN_term_num=(dx_dt_1PN_term_base.numerator().collect(e_t))*x^6
dx_dt_1PN_term_den=dx_dt_1PN_term_base.denominator()
dx_dt_1PN_term=dx_dt_1PN_term_num/dx_dt_1PN_term_den


# In[23]:


dx_dt_15PN_term_base=dx_dt.coefficient(Xi,3).coefficient(x^(13/2)).subs({e_t^10 - 5*e_t^8 + 10*e_t^6 - 10*e_t^4 + 5*e_t^2 - 1:-(-e_t^2+1)^5})
dx_dt_15PN_term_num=(dx_dt_15PN_term_base.numerator().collect(e_t))*x^(13/2)
dx_dt_15PN_term_den=dx_dt_15PN_term_base.denominator()
dx_dt_15PN_term=dx_dt_15PN_term_num/dx_dt_15PN_term_den


# In[24]:


dx_dt_2PN_term_base_0=(dx_nonspinning_dt.coefficient(Xi,4)+dx_dt.coefficient(Xi,4).coefficient(SO,1)*SO+dx_dt.coefficient(Xi,4).coefficient(SO,2)*SO^2)
dx_dt_2PN_term_base=(dx_dt_2PN_term_base_0.coefficient(x^7).subs({(-e_t^2+1):y}).subs({(e_t^2-1):-y}))
dx_dt_2PN_term_pt_1=(dx_dt_2PN_term_base.coefficient(y^(-5)).numerator().collect(e_t)/dx_dt_2PN_term_base.coefficient(y^(-5)).denominator())*(1/y^5)
dx_dt_2PN_term_pt_2=(dx_dt_2PN_term_base.coefficient(y^(-11/2)).numerator().collect(e_t)/dx_dt_2PN_term_base.coefficient(y^(-11/2)).denominator())*(1/y^(11/2))
dx_dt_2PN_term=((dx_dt_2PN_term_pt_1+dx_dt_2PN_term_pt_2).subs({y:-e_t^2+1}))*x^7


# In[25]:


dx_dt_25PN_term_base_0=(dx_dt.coefficient(Xi,5)+dx_dt_horiz_labelled.coefficient(Xi,5)).coefficient(x^(15/2)).subs({(-e_t^2+1)^(1/2):y^(1/2)}).subs({e_t^12 - 6*e_t^10 + 15*e_t^8 - 20*e_t^6 + 15*e_t^4 - 6*e_t^2 + 1:y^6})
dx_dt_25PN_term_pt_1=(dx_dt_25PN_term_base_0.coefficient(y^(-6)).numerator().collect(e_t)/dx_dt_25PN_term_base_0.coefficient(y^(-6)).denominator())*(1/y^6)
dx_dt_25PN_term_pt_2=(dx_dt_25PN_term_base_0.coefficient(y^(-11/2)).numerator().collect(e_t)/dx_dt_25PN_term_base_0.coefficient(y^(-11/2)).denominator())*(1/y^(11/2))


# In[26]:


dx_dt_25PN_term_base_1=(dx_nonspinning_dt.coefficient(Xi,5).coefficient(x^(15/2)).subs({-e_t^2+1:y}).subs({e_t^2-1:-y}))
dx_dt_25PN_term_pt_3=(dx_dt_25PN_term_base_1.subs({y:1/z}).coefficient(z,0).numerator().collect(e_t)/dx_dt_25PN_term_base_1.subs({y:1/z}).coefficient(z,0).denominator())
dx_dt_25PN_term_pt_4=(dx_dt_25PN_term_base_1.coefficient(y^(-1)).numerator().collect(e_t)/dx_dt_25PN_term_base_1.coefficient(y^(-1)).denominator())*(1/y)
dx_dt_25PN_term_pt_5=(dx_dt_25PN_term_base_1.coefficient(y^(-3/2)).numerator().collect(e_t)/dx_dt_25PN_term_base_1.coefficient(y^(-3/2)).denominator())*(1/y^(3/2))
dx_dt_25PN_term=((dx_dt_25PN_term_pt_1+dx_dt_25PN_term_pt_2+dx_dt_25PN_term_pt_3+dx_dt_25PN_term_pt_4+dx_dt_25PN_term_pt_5).subs({y:-e_t^2+1}))*x^(15/2)


# In[27]:


dx_dt_3PN_term_base_0=((dx_dt.coefficient(Xi,6)+dx_nonspinning_dt.coefficient(Xi,6)).expand().coefficient(x^8).subs({e_t^12 - 6*e_t^10 + 15*e_t^8 - 20*e_t^6 + 15*e_t^4 - 6*e_t^2 + 1:y^6}).subs({(-e_t^2+1):y}))
dx_dt_3PN_term_pt_1=(dx_dt_3PN_term_base_0.coefficient(y^(-6)).numerator().collect(e_t)/dx_dt_3PN_term_base_0.coefficient(y^(-6)).denominator())*(1/y^6)
dx_dt_3PN_term_pt_2=(dx_dt_3PN_term_base_0.coefficient(y^(-13/2)).numerator().collect(e_t)/dx_dt_3PN_term_base_0.coefficient(y^(-13/2)).denominator())*(1/y^(13/2))
dx_dt_3PN_term_pt_3=(dx_dt_3PN_term_base_0.subs({y:1/z}).coefficient(z,0).numerator().collect(e_t)/dx_dt_3PN_term_base_0.subs({y:1/z}).coefficient(z,0).denominator())
dx_dt_3PN_term=((dx_dt_3PN_term_pt_1+dx_dt_3PN_term_pt_2+dx_dt_3PN_term_pt_3).subs({y:-e_t^2+1}))*x^8


# In[28]:


dx_dt_final=dx_dt_N_term+dx_dt_1PN_term*Xi^2+dx_dt_15PN_term*Xi^3+dx_dt_2PN_term*Xi^4+dx_dt_25PN_term*Xi^5+dx_dt_3PN_term*Xi^6


# **Factor the de_dt expressions and include non spinning terms from Michael, Faye,...**

# In[29]:


de_dt_N_term=(de_dt.coefficient(Xi,0).subs({(-e_t^2+1)^(-5/2):y}).factor().subs({y:(-e_t^2+1)^(-5/2)}))


# In[30]:


de_dt_1PN_term_base=(de_dt.coefficient(Xi,2).coefficient(x^5).subs({(-e_t^2+1)^(-7/2):y}).factor().subs({y:(-e_t^2+1)^(-7/2)}))
de_dt_1PN_term_num=(de_dt_1PN_term_base.numerator().collect(e_t))*x^5
de_dt_1PN_term_den=de_dt_1PN_term_base.denominator()
de_dt_1PN_term=de_dt_1PN_term_num/de_dt_1PN_term_den


# In[31]:


de_dt_15PN_term_base=de_dt.coefficient(Xi,3).coefficient(x^(11/2)).subs({e_t^8 - 4*e_t^6 + 6*e_t^4 - 4*e_t^2 + 1:(-e_t^2+1)^4})
de_dt_15PN_term_num=(de_dt_15PN_term_base.numerator().collect(e_t))*x^(11/2)
de_dt_15PN_term_den=de_dt_15PN_term_base.denominator()
de_dt_15PN_term=de_dt_15PN_term_num/de_dt_15PN_term_den


# In[32]:


de_dt_2PN_term_base=de_dt.coefficient(Xi,4).coefficient(x^6).subs({e_t^8 - 4*e_t^6 + 6*e_t^4 - 4*e_t^2 + 1:y^4}).subs({(-e_t^2+1):y})
de_dt_2PN_term_pt_1=(de_dt_2PN_term_base.coefficient(y^(-9/2)).numerator().collect(e_t)/de_dt_2PN_term_base.coefficient(y^(-9/2)).denominator())*(1/y^(9/2))
de_dt_2PN_term_pt_2=(de_dt_2PN_term_base.coefficient(y^(-4)).numerator().collect(e_t)/de_dt_2PN_term_base.coefficient(y^(-4)).denominator())*(1/y^(4))
de_dt_2PN_term=((de_dt_2PN_term_pt_1+de_dt_2PN_term_pt_2).subs({y:-e_t^2+1}))*x^6


# In[33]:


de_dt_25PN_term_base_0=(de_dt.coefficient(Xi,5)+de_dt_horiz_labelled.coefficient(Xi,5)).coefficient(x^(13/2)).subs({e_t^10 - 5*e_t^8 + 10*e_t^6 - 10*e_t^4 + 5*e_t^2 - 1:-y^5}).subs({(-e_t^2+1):y})
de_dt_25PN_term_pt_1=(de_dt_25PN_term_base_0.coefficient(y^(-5)).numerator().collect(e_t)/de_dt_25PN_term_base_0.coefficient(y^(-5)).denominator())*(1/y^5)
de_dt_25PN_term_pt_2=(de_dt_25PN_term_base_0.coefficient(y^(-9/2)).numerator().collect(e_t)/de_dt_25PN_term_base_0.coefficient(y^(-9/2)).denominator())*(1/y^(9/2))


# In[34]:


de_dt_25PN_term_base_1=(de_nonspinning_dt.coefficient(Xi,5).coefficient(x^(13/2)).subs({-e_t^2+1:y}).subs({e_t^2-1:-y}))
de_dt_25PN_term_pt_3=(de_dt_25PN_term_base_1.coefficient(y^(-1/2)).numerator().collect(e_t)/de_dt_25PN_term_base_1.coefficient(y^(-1/2)).denominator())*(1/y^(1/2))
de_dt_25PN_term_pt_4=(de_dt_25PN_term_base_1.coefficient(y^(1/2)).numerator().collect(e_t)/de_dt_25PN_term_base_1.coefficient(y^(1/2)).denominator())*(y^(1/2))
de_dt_25PN_term_pt_5=(de_dt_25PN_term_base_1.coefficient(y,0).numerator().collect(e_t)/de_dt_25PN_term_base_1.coefficient(y,0).denominator())
de_dt_25PN_term=((de_dt_25PN_term_pt_1+de_dt_25PN_term_pt_2+de_dt_25PN_term_pt_3+de_dt_25PN_term_pt_4+de_dt_25PN_term_pt_5).subs({y:-e_t^2+1}))*x^(13/2)


# In[35]:


de_dt_3PN_term_base_0=((de_dt.coefficient(Xi,6)+de_nonspinning_dt.coefficient(Xi,6)).expand().coefficient(x^7).subs({e_t^10 - 5*e_t^8 + 10*e_t^6 - 10*e_t^4 + 5*e_t^2 - 1:-y^5}).subs({(-e_t^2+1):y}).subs({(e_t^2-1):-y}))
de_dt_3PN_term_base_1=(de_dt_3PN_term_base_0-(de_dt_3PN_term_base_0.coefficient((sqrt(y)+1)^(-1)))*((sqrt(y)+1)^(-1))).expand()
de_dt_3PN_term_base_2=de_dt_3PN_term_base_1+(de_dt_3PN_term_base_0.coefficient((sqrt(y)+1)^(-1)).subs({y^(-11/2):y_0}))*((sqrt(y)+1)^(-1))


# In[36]:


de_dt_3PN_term_pt_1=(de_dt_3PN_term_base_2.coefficient(y^(-5)).numerator().collect(e_t)/de_dt_3PN_term_base_2.coefficient(y^(-5)).denominator())*(1/y^5)
de_dt_3PN_term_pt_2=((de_dt_3PN_term_base_2.coefficient(y^(-11/2)).numerator().collect(e_t)/de_dt_3PN_term_base_2.coefficient(y^(-11/2)).denominator())+de_dt_3PN_term_base_2.coefficient(y_0))*(1/y^(11/2))
de_dt_3PN_term_pt_3=(de_dt_3PN_term_base_2.coefficient(y^(1/2)).numerator().collect(e_t)/de_dt_3PN_term_base_2.coefficient(y^(1/2)).denominator())*(y^(1/2))


# In[37]:


messy_term=(de_dt_3PN_term_base_2.coefficient(y,0)-(de_dt_3PN_term_base_2.coefficient(y_0))*(y_0)).expand()


# In[38]:


de_dt_3PN_term_pt_4=(messy_term.coefficient(y,0).numerator().collect(e_t)/messy_term.coefficient(y,0).denominator())
de_dt_3PN_term=((de_dt_3PN_term_pt_1+de_dt_3PN_term_pt_2+de_dt_3PN_term_pt_3+de_dt_3PN_term_pt_4).subs({y:-e_t^2+1}))*x^7


# In[39]:


de_dt_final=de_dt_N_term+de_dt_1PN_term*Xi^2+de_dt_15PN_term*Xi^3+de_dt_2PN_term*Xi^4+de_dt_25PN_term*Xi^5+de_dt_3PN_term*Xi^6


# **Newtonian**

# In[40]:


dx_de_N_integrand=(dx_dt_N_term/de_dt_N_term)/x


# In[41]:


dx_de_N=dx_de_N_integrand.integrate(e_t).subs({- 6/19*log(e_t^2):log(e_t^(-12/19))}).subs({-870/2299*log(121*e_t^2 + 304):log((121*e_t^2 + 304)^(-870/2299))})


# In[42]:


dx_de_N_2=dx_de_N-dx_de_N.subs({e_t:e0})


# In[43]:


x_N=(exp(dx_de_N_2).simplify_full())


# **1 PN**

# In[44]:


dx_de_1PN_integrand=(((dx_dt_N_term+dx_dt_1PN_term*Xi^2)/(de_dt_N_term+de_dt_1PN_term*Xi^2))/x)


# For the step below you can also just do dx_de_1PN_integrand.taylor(Xi,0,2) and then do .coefficient(Xi,2).taylor(e_t,...))

# In[45]:


dx_de_1PN_integrand_2=dx_de_1PN_integrand.taylor((Xi,0),(e_t,0),8).coefficient(Xi,2)


# In[46]:


dx_de_1PN_integrand_3=dx_de_1PN_integrand_2.factor().subs({x:x_N_1}).taylor((e_t,0),(e0,0),6)


# In[47]:


dx_de_1PN=(dx_de_1PN_integrand_3.integrate(e_t))*Xi^2


# In[48]:


dx_de_1PN_2=dx_de_1PN-dx_de_1PN.subs({e_t:e0})


# In[49]:


x_1PN=((exp(dx_de_1PN_2).simplify_full())*x_N).simplify_full()


# **1.5 PN**

# In[50]:


dx_de_15PN_integrand=(((dx_dt_N_term+dx_dt_1PN_term*Xi^2+dx_dt_15PN_term*Xi^3)/(de_dt_N_term+de_dt_1PN_term*Xi^2+de_dt_15PN_term*Xi^3))/x)


# In[51]:


dx_de_15PN_integrand_2=dx_de_15PN_integrand.taylor((Xi,0),(e_t,0),9).coefficient(Xi,3)


# In[52]:


dx_de_15PN_integrand_3=(dx_de_15PN_integrand_2.collect(x).subs({x:x_N_1})).taylor((e_t,0),(e0,0),6)


# In[53]:


dx_de_15PN=(dx_de_15PN_integrand_3.integrate(e_t))*Xi^3


# In[54]:


dx_de_15PN_2=dx_de_15PN-dx_de_15PN.subs({e_t:e0})


# In[55]:


x_15PN=(((exp(dx_de_15PN_2).simplify_full())*x_1PN).simplify_full())


# **2 PN**

# In[56]:


dx_de_2PN_integrand=(((dx_dt_N_term+dx_dt_1PN_term*Xi^2+dx_dt_15PN_term*Xi^3+dx_dt_2PN_term*Xi^4)/(de_dt_N_term+de_dt_1PN_term*Xi^2+de_dt_15PN_term*Xi^3+de_dt_2PN_term*Xi^4))/x)


# In[57]:


dx_de_2PN_integrand_2=dx_de_2PN_integrand.taylor((Xi,0),(e_t,0),10).coefficient(Xi,4)


# In[58]:


dx_de_2PN_integrand_3=dx_de_2PN_integrand_2.collect(x).subs({x:x_N_1}).taylor((e_t,0),(e0,0),6)


# In[59]:


dx_de_2PN_integrand_4=(dx_de_1PN_integrand_2.collect(x).subs({x:x_N_1+Xi^2*x_1PN_final}).taylor((Xi,0),(e_t,0),(e0,0),10).coefficient(Xi,2))*Xi^4


# In[60]:


dx_de_2PN=(dx_de_2PN_integrand_3.integrate(e_t)+dx_de_2PN_integrand_4.coefficient(Xi,4).expand().integrate(e_t))*Xi^4


# In[61]:


dx_de_2PN_2=dx_de_2PN-dx_de_2PN.subs({e_t:e0})


# In[62]:


x_2PN=(((exp(dx_de_2PN_2).simplify_full())*x_15PN).simplify_full())


# **2.5 PN**

# In[63]:


dx_de_25PN_integrand=(((dx_dt_N_term+dx_dt_1PN_term*Xi^2+dx_dt_15PN_term*Xi^3+dx_dt_2PN_term*Xi^4+dx_dt_25PN_term*Xi^5)/(de_dt_N_term+de_dt_1PN_term*Xi^2+de_dt_15PN_term*Xi^3+de_dt_2PN_term*Xi^4+de_dt_25PN_term*Xi^5))/x)


# In[64]:


dx_de_25PN_integrand_2=dx_de_25PN_integrand.taylor((Xi,0),(e_t,0),11).coefficient(Xi,5)


# In[65]:


dx_de_25PN_integrand_3=dx_de_25PN_integrand_2.collect(x).subs({x:x_N_1}).taylor((e_t,0),(e0,0),6)


# In[66]:


dx_de_25PN_integrand_4=dx_de_1_contribution.coefficient(Xi,5)


# In[67]:


dx_de_25PN_integrand_5=dx_de_15_contribution.coefficient(Xi,5)


# In[68]:


dx_de_25PN_pt1=dx_de_25PN_integrand_3.integrate(e_t)


# In[69]:


dx_de_25PN_pt2=dx_de_25PN_integrand_4.expand().integrate(e_t)


# In[70]:


dx_de_25PN_pt3=dx_de_25PN_integrand_5.expand().integrate(e_t)


# In[71]:


dx_de_25PN=(dx_de_25PN_pt1+dx_de_25PN_pt2+dx_de_25PN_pt3)*Xi^5


# In[72]:


dx_de_25PN_2=dx_de_25PN-dx_de_25PN.subs({e_t:e0})


# In[73]:


x_25PN=(((exp(dx_de_25PN_2).simplify_full())*x_2PN).simplify_full())


# In[74]:


x_25PN_1=(x_25PN.taylor((Xi,0),(e_t,0),(e0,0),13).coefficient(Xi,5))


# In[75]:


x_25PN_final=((x_25PN_1*(e_factor)).expand().taylor((e_t,0),(e0,0),6))*(x0/e_factor)


# In[76]:


filename = "../x_e_files/x_25PN.txt"


with open(filename, "w") as f:
    f.write("# Expansion of x_25PN\n")
    f.write("# Plain ASCII output from Sage\n\n")
    f.write("x_25PN =\n")
    f.write(str(x_2PN_expr+x_25PN_final*Xi^5))
    f.write("\n\n")


# In[80]:


#Labelling this x_2PN caused huge error. Be careful. Renamed to x_2PN_expr
filename = "../x_e_files/x_25PN.txt"

with open(filename, "r") as f:
    content = f.read()


match = re.search(r"x_25PN\s*=\s*(.*)", content, re.DOTALL)
if match:
    expr_str = match.group(1).strip()
    # Convert string back into a Sage expression
    x_25PN_expr = (sage_eval(expr_str, locals=globals())).expand()


# **3 PN**

# In[75]:


dx_de_3PN_integrand=(((dx_dt_N_term+dx_dt_1PN_term*Xi^2+dx_dt_15PN_term*Xi^3+dx_dt_2PN_term*Xi^4+dx_dt_25PN_term*Xi^5+dx_dt_3PN_term*Xi^6)/(de_dt_N_term+de_dt_1PN_term*Xi^2+de_dt_15PN_term*Xi^3+de_dt_2PN_term*Xi^4+de_dt_25PN_term*Xi^5+de_dt_3PN_term*Xi^6))/x)


# In[76]:


dx_de_3PN_integrand_2=dx_de_3PN_integrand.taylor((Xi,0),(e_t,0),12).coefficient(Xi,6)


# In[77]:


dx_de_3PN_integrand_3=dx_de_3PN_integrand_2.collect(x).subs({x:x_N_1}).taylor((e_t,0),(e0,0),6)


# In[78]:


dx_de_3PN_integrand_4=dx_de_1_contribution.coefficient(Xi,6)


# In[79]:


dx_de_3PN_integrand_5=dx_de_15_contribution.coefficient(Xi,6)


# In[80]:


dx_de_3PN_integrand_6=dx_de_2_contribution.coefficient(Xi,6)


# In[81]:


dx_de_3PN_pt1=dx_de_3PN_integrand_3.integrate(e_t)


# In[82]:


dx_de_3PN_pt2=dx_de_3PN_integrand_4.expand().integrate(e_t)


# In[83]:


dx_de_3PN_pt3=dx_de_3PN_integrand_5.expand().integrate(e_t)


# In[84]:


dx_de_3PN_pt4=dx_de_3PN_integrand_6.expand().integrate(e_t)


# In[85]:


dx_de_3PN=(dx_de_3PN_pt1+dx_de_3PN_pt2+dx_de_3PN_pt3+dx_de_3PN_pt4)*Xi^6


# In[86]:


dx_de_3PN_2=dx_de_3PN-dx_de_3PN.subs({e_t:e0})


# In[ ]:


x_3PN=(((exp(dx_de_3PN_2).simplify_full())*x_25PN).simplify_full())


# In[ ]:


x_3PN_1=(x_3PN.taylor((Xi,0),(e_t,0),(e0,0),14).coefficient(Xi,6))


# In[ ]:


x_3PN_final=((x_3PN_1*(e_factor)).expand().taylor((e_t,0),(e0,0),6))*(x0/e_factor)


# In[ ]:


filename = "../x_e_files/x_3PN.txt"


with open(filename, "w") as f:
    f.write("# Expansion of x_3PN\n")
    f.write("# Plain ASCII output from Sage\n\n")
    f.write("x_3PN =\n")
    f.write(str(x_25PN_expr+x_3PN_final*Xi^6))
    f.write("\n\n")


# In[87]:


#show(((x_25PN_expr*e_factor).expand().coefficient(Xi,5).expand().coefficient(SO,3)-xe_expr.coefficient(Xi,5).coefficient(SO,3)).expand())

