# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 09:21:55 2019

@author: Ricardo Luna Hernández

pyporap: Post-Regession Analysis Package.

This package contains routines for post-regression analysis:
sensitivity, identifiability, and significance of parameters
in a dynamic model. This package is based on the CasADi package,
an open-source optimization software that uses automatic differentiation.
"""
# %% importar librerías
from casadi import*
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.integrate import simpson
# %% Analisis de sensibilidad
def sensitivity_analysis(x,p,ode,x0,tspan,param,solver,stype):
    """ This function performs sensitivity analysis of a dynamic model 
    and returns sensitivity plots over time and the accumulated
    sensitivity integral.
    Inputs:
        x: vector de estados del modelo (Casadi symbolic)
        p: vector de parámetros del modelo (Casadi symbolic)
        ode: vector de ecuaciones diferenciales del modelo (Casadi symbolic)
        x0: vector de condiciones iniciales del modelo (list)
        tspan: vector de tiempo de integración (array)
        param: vector de parámetros del modelo a evaluar (list)
        solver: solver de integración (str: 'idas', 'cvodes')
        stype: tipo de sensibilidad (str: 'abs', 'absp', 'rel1', 'rel2')
                'abs': sensisibilidad absoluta
                'absp': sensibilidad * valor nomimal parámetro
                'rel1': sensibilidad * (valor nominal parámetro/valor estado)
                'rel2': sensibilidad * (valor nominal parámetro/max valor estado)
                'norm': valor absoluto |sensibilidad * (valor nominal parámetro/max valor estado)|
                'mean': sensibilidad normalizada relativa promedio en el tiempo
        """    
    n_state  = x.shape[0] # número de estados
    n_par    = p.shape[0] # número de parámetros
    G,xstate = sensitivity_ode_integration(x,p,ode,x0,tspan,param,solver)        
    G_mean          = np.zeros([n_state,n_par]) # matriz de sensibilidad
    names           = [0]*n_par # etiqueta del número de parámetro    
    for i in np.arange(0,n_state):       
        for j in np.arange(0,n_par):
            if stype == 'abs':
                plt.figure(i)
                plt.subplot(math.ceil(n_par/2),2,j+1)
                plt.plot(tspan,G[j*n_state+i,:])
                plt.xlabel('Time',fontsize=14)
                plt.ylabel(r'$dx_{'+str(i+1)+'}$'+r'$/dp_{'+str(j+1)+'}$',\
                                  fontsize=14)
            elif stype == 'absp':
                plt.figure(i)
                plt.subplot(math.ceil(n_par/2),2,j+1)
                plt.plot(tspan,G[j*n_state+i,:]*(param[j]))
                plt.xlabel('Time',fontsize=14)
                plt.ylabel(r'$p_{'+str(j+1)+'}$'+r'$\cdot dx_{'+str(i+1)+'}$'+\
                                 r'$/dp_{'+str(j+1)+'}$',fontsize=14)
            elif stype == 'rel1':
                plt.figure(i)
                plt.subplot(math.ceil(n_par/2),2,j+1)
                plt.plot(tspan,G[j*n_state+i,:]*(param[j]/xstate[i,:]))
                plt.xlabel('Time',fontsize=14)
                plt.ylabel(r'$(p_{'+str(j+1)+'}$'+r'$/x_{'+str(i+1)+'})$'+\
                           r'$\cdot dx_{'+str(i+1)+'}$'+r'$/dp_{'+str(j+1)+'}$',\
                           fontsize=14)
            elif stype == 'rel2':
                plt.figure(i)
                plt.subplot(math.ceil(n_par/2),2,j+1)
                plt.plot(tspan,G[j*n_state+i,:]*(param[j]/np.max(xstate[i,:])))
                plt.xlabel('Time',fontsize=14)
                plt.ylabel(r'$(p_{'+str(j+1)+'}$'+r'$/max(x_{'+str(i+1)+'}))$'+\
                           r'$\cdot dx_{'+str(i+1)+'}$'+r'$/dp_{'+str(j+1)+'}$',\
                           fontsize=11)
            elif stype == 'norm':
                plt.figure(i)
                plt.subplot(math.ceil(n_par/2),2,j+1)
                plt.plot(tspan,np.abs(G[j*n_state+i,:]*(param[j]/np.max(xstate[i,:]))))
                plt.xlabel('Time',fontsize=14)
                plt.ylabel(r'$|$'+r'$(p_{'+str(j+1)+'}$'+r'$/max(x_{'+str(i+1)+'}))$'+\
                           r'$\cdot dx_{'+str(i+1)+'}$'+r'$/dp_{'+str(j+1)+'}$'+r'$|$',\
                           fontsize=11)
            elif stype == 'mean':
                plt.figure(n_state)
            G_mean[i,j] = simpson(np.abs(G[j*n_state+i,:]*\
                  (param[j]/np.max(xstate[i,:])))*(1/tspan[-1]),tspan)
            names[j]=r'$p_{'+str(j+1)+'}$'       
        plt.tight_layout(h_pad=0.5)
        if stype == 'mean':
            plt.figure(n_state)
            plt.subplot(n_state,1,i+1)
            plt.bar(names,G_mean[i,:])
            plt.ylabel(r'$x_{'+str(i+1)+'}$')
    if stype == 'mean':
        return G_mean, None
    
    return G, xstate
# %% Ecuaciones diferenciales de sensibilidad
def sensitivity_ode(x,p,ode):
    """ Esta función crea el ode de las ecuaciones de sensibilidad parámetrica 
    acopladas al ode del modelo.
    inputs:
        x: vector de estados del modelo (Casadi symbolic)
        p: vector de parámetros del modelo (Casadi symbolic)
        ode: vector de ecuaciones diferenciales del modelo (Casadi symbolic)"""    
    n_state     = x.shape[0] # número de estados
    n_par       = p.shape[0] # número de parámetros
    G           = MX.sym('G',n_state,n_par) # se crea variable G (matriz de sensibilidad)
    xG          = vertcat(x,reshape(G,n_state*n_par,1)) # concatenación de estados

    # del modelo con estados nuevos estados asociados a la sensibilidad    
    dfdp        = jacobian(ode,p) # calculo del jacobiano df/dp
    dfdx        = jacobian(ode,x) # calculo del jacobiano df/dx
    dGdt        = (dfdx@G)+dfdp   # ecuaciones diferenciales de la sensibilidad
    dGdt        = reshape(dGdt,n_state*n_par,1) 
    xGdot       = vertcat(ode,dGdt) # concatenación de ecuaciones 

    # diferenciales del modelo con ecuaciones difereciales de sensibilidad    
    return xG, xGdot;
# %% Integración en el tiempo del modelo + sensibilidad
def sensitivity_ode_integration(x,p,ode,x0,tspan,param,solver):
    """ Esta función realiza la integración del modelo acoplado al nuevo ode 
    de las ecuaciones diferenciales de sensibilidad.
    Inputs:
        x: vector de estados del modelo (Casadi symbolic)
        p: vector de parámetros del modelo (Casadi symbolic)
        ode: vector de ecuaciones diferenciales del modelo (Casadi symbolic)
        x0: vector de condiciones iniciales del modelo (list)
        tspan: vector de tiempo de integración (array)
        param: vector de parámetros del modelo a evaluar (list) 
        solver:solver de integración (str: 'idas', 'cvodes')"""
    xG, xGdot       = sensitivity_ode(x,p,ode) # se llama a la función
    # sensitivity_ode para acoplar el ode del modelo con ode de sensibilidad
    dae         = {'x':xG,'ode':xGdot,'z':[],'alg':[],'p':p} #creación del dae
    opts        = {'tf':tspan[-1],'grid':tspan,'output_t0':True} # opciones de integración
    Sen_time    = integrator('Sensitivity_time',solver,dae,opts) # se crea el objeto de integración
    n_state     = x.shape[0] # número de estados
    n_par       = p.shape[0] # número de parámetros
    G           = np.zeros([1,n_state*n_par]) # se crea G (matriz de sensibilidad)
    x0_new      = x0+[0]*n_state*n_par; # condiciones iniciales del ode
    result      = Sen_time(x0=x0_new,p=param) # se ejecuta el integrador
    S           = result['xf'].full() # se extraen los resultados
    xstate      = S[0:n_state,:] # variables de estado del modelo
    G           = S[n_state:,:]  # variables de estado de la sensibilidad
    return G,xstate
# %% Análisis de identificabilidad    
def identifiability(x,p,ode,x0,tspan,param,solver,Cmax):
    """ 
    Esta función realiza el análisis de identificabilidad de los parámetros del
    modelo dinámico y entrega como resultados la matriz de correlación de los
    parámetros.
    Inputs:
        x: vector de estados del modelo (Casadi symbolic)
        p: vector de parámetros del modelo (Casadi symbolic)
        ode: vector de ecuaciones diferenciales del modelo (Casadi symbolic)
        x0: vector de condiciones iniciales del modelo (list)
        tspan: vector de tiempo de integración (array)
        param: vector de parámetros del modelo a evaluar (list)
        solver:solver de integración (str: 'idas', 'cvodes')
        Cmax: umbral de correlación de los parámetros (0.95 valor tipico usado)
    """
    n_state = x.shape[0] # número de estados
    n_par   = p.shape[0] # número de parámetros
    n_time  = len(tspan) # pasos de tiempo
    G,xstate = sensitivity_ode_integration(x,p,ode,x0,tspan,param,solver) 
    Gaux            = reshape(G.T,n_time*n_state, n_par)
    C               = np.corrcoef(Gaux.T)
    plt.figure()    
    plt.pcolor(C,cmap='seismic')
    plt.colorbar()    
    names = [0]*n_par    
    for i in np.arange(n_par):
        names[i]=r'$p_{'+str(i+1)+'}$'
    plt.xticks(np.arange(n_par)+0.5,names,fontsize=12)
    plt.yticks(np.arange(n_par)+0.5,names,fontsize=12)
    plt.show()
    print('parameters with high correlation >=',Cmax)
    for i in np.arange(n_par):
        for j in np.arange(n_par):
            if np.abs(C[i,j]) >= Cmax and i!=j:
                print(i+1,"and",j+1)
    return C 
# %% Análisis de significancia  
def significance(x,p,ode,x0,tspan,param,solver,sigma2_x):
    """ 
    Esta función realiza el análisis de significancia de los parámetros del
    modelo dinámico. 
    Inputs:
        x: vector de estados del modelo (Casadi symbolic)
        p: vector de parámetros del modelo (Casadi symbolic)
        ode: vector de ecuaciones diferenciales del modelo (Casadi symbolic)
        x0: vector de condiciones iniciales del modelo (list)
        tspan: vector de tiempo de integración (array)
        param: vector de parámetros del modelo a evaluar (list)
        solver:solver de integración (str: 'idas', 'cvodes')
        sigma2_x: vector varianza de la medición de los estados para construir
        matrix de covarianza.
    Outputs:
        sigma: vector de desviación estándar de los parámetros
        t_value: student t-value (parametro/desviación estándar parámetro)
        CI_l: intervalo de confianza límite inferior (95%) 
        CI_u: intervalo de confianza límite superior (95%)
    """
    n_state = x.shape[0] # número de estados
    n_par   = p.shape[0] # número de parámetros
    n_time  = len(tspan) # pasos de tiempo    
    G,xstate = sensitivity_ode_integration(x,p,ode,x0,tspan,param,solver)
    FIM = np.zeros([n_par,n_par]) # matriz de fisher
    Q   = np.eye(n_state)*(1/sigma2_x) # Inverso matriz covarianza
    G_old = np.zeros([n_state,n_par]) # matriz de sensibilidad 
    
    for k in np.arange(0,n_time):        
        for i in np.arange(0,n_state):
            for j in np.arange(0,n_par):
                G_old[i,j]= G[i+j*n_state,k] # se rellena matriz de sensibilidad       
        FIM = G_old.T@Q@G_old+FIM # se calcula la matriz de Fisher
    sigma2   = np.diagonal(np.linalg.inv(FIM)) # se estima la varianza de los 
    # parámetros como el inverso de la matriz de fisher (elementos de la diagonal)
    sigma    = np.sqrt(sigma2) # se calcula la desviación estándar de los parámetros
    t_value  = param/np.sqrt(sigma2) # se calcula el t-value (sacher et al. 2011)
    CI_l     = param-2*sigma # intervalo de confianza inferior (95%)
    CI_u     = param+2*sigma # intervalo de confianza superior (95%)
    #CC       = 3.92*(sigma/param)  # sanchez
    #CV       = (sigma/param)*100    
    return sigma, t_value, CI_l, CI_u