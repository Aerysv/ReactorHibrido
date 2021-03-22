/*-----------------------------------------------------------------------------------------
 LIBRARY: AJUSTE_PARAMETROS
 FILE: Ajuste_8Feb_2021
 CREATION DATE: 21/03/2021
 DESCRIPTION:
 Estados medidos:
 	- Ca, Cb, Cc, Cd
	- T
	- Tc
Variables de decisión:
	- Ca_ini, Cb_ini, Cc_ini, Cd_ini
	- T_ini
	- Tc_ini
	- k10, k20, k30
	- Ea1, Ea2, Ea3
	- dHrxn1, dHrxn2, dHrxn3
	
-----------------------------------------------------------------------------------------*/
USE MATH
COMPONENT Ajuste_6estados (INTEGER nsal = 6)

	DATA

		REAL V = 11.5			"Volumen del reactor (L)"
		REAL Vc = 1				"Volumen de la camisa (L)"
		REAL R = 0.00831		"Constante de los gases (kJ/mol/K)"
		REAL rho = 1			"Densidad del agua (kg/L)"
		REAL Cp = 4.184		"Capacidad Calorifica del agua (kJ/kg/C)"	
      REAL Ca0 = 5			"Concentracion de entrada de A (mol/L)"
		
   	-- Variables de decisión
		REAL k10_estim = 6.4e9				"Factor preexponencial (1/min)"
		REAL k20_estim = 4.84e10			"Factor preexponencial (1/min)"
		REAL k30_estim = 8.86e11			"Factor preexponencial (1/min/mol)"
		REAL Ea1 = 52.1				"Energía de activacion kJ/mol"
		REAL Ea2 = 70.0				"Energía de activacion kJ/mol"
		REAL Ea3 = 65.0				"Energía de activacion kJ/mol"
		REAL dHrxn1 = -28.75			"Calor de reaccion (kJ/mol)"
		REAL dHrxn2 = -34.50			"Calor de reaccion (kJ/mol)"
		REAL dHrxn3 = -40.25			"Calor de reaccion (kJ/mol)"
		REAL alpha = 1.2197
		REAL y_ini[6] = {23.5, 14, 3.5, 0.4, 0.11, 0.57}
	
		-- Optimizador
		REAL pesos[6] = {100, 100, 100, 100, 100, 100}
		REAL media[6] = {0.104, 3.44, 0.78, 0.34, 24.9, 12.3}
		REAL LiminfT = 5
    	REAL LimsupT = 60
    	REAL LiminfTc = 5
    	REAL LimsupTc = 60
	 	REAL LiminfCa = 0
    	REAL LimsupCa = 5
		REAL C = 5
		
	DECLS
		-- Caudales
		REAL q                "Caudal reactivos (L/min)"
		REAL Fr               "Caudal refrigerante (L/min)"
	
    	-- Salidas medidas del proceso
		REAL T                 "Temperatura reactor (C)"
		REAL T0					  "Temperatura entrada reactivos (C)"
		REAL Tc                "Temperatura refrigerante (C)"
		REAL Tc0					  "Temperatura entrada refrigerante (C)"

		-- Otras variables del modelo
		REAL Q					"Calor transferido (kJ/mol)"
		REAL Heat_rxn			"Calor generado (kJ/min)"
		REAL Ca, Cb, Cc, Cd
		REAL r1, r2, r3
		REAL UA
		
		-- Lectura de datos
		TABLE_1D tab_q, tab_Fr, tab_T, tab_Tc, tab_T0, tab_Tc0, tab_Ca, tab_Cb, tab_Cc, tab_Cd
		
		-- Optimizador
		REAL F_optim[7]
		REAL J[nsal]
		REAL J_costo
		REAL y_real[nsal], y_modelo[nsal]
		REAL k10, k20, k30
		
	INIT

		-- Lectura de variables
		-- Variables de entrada
		q = linearInterp1D(tab_q, TIME)
		Fr = linearInterp1D(tab_Fr, TIME)
		
		-- Perturbaciones
		T0 = linearInterp1D(tab_T0, TIME)
		Tc0 = linearInterp1D(tab_Tc0, TIME)

		-- Variables de salida
		y_real[1] = linearInterp1D(tab_Ca, TIME)
		y_real[2] = linearInterp1D(tab_Cb, TIME)
		y_real[3] = linearInterp1D(tab_Cc, TIME)
		y_real[4] = linearInterp1D(tab_Cd, TIME)
		y_real[5] = linearInterp1D(tab_T, TIME)
		y_real[6] = linearInterp1D(tab_Tc, TIME)

		y_modelo[1] = y_real[1]
		y_modelo[2] = y_real[2]
		y_modelo[3] = y_real[3]
		y_modelo[4] = y_real[4]
		y_modelo[5] = y_real[5]
		y_modelo[6] = y_real[6]


		Ca = y_modelo[1]
		Cb = y_modelo[2]
		Cc = y_modelo[3]
		Cd = y_modelo[4]
		T  = y_modelo[5]
		Tc = y_modelo[6]		
		/*
		
		y_modelo[1] = y_ini[1]
		y_modelo[2] = y_ini[2]
		y_modelo[3] = y_ini[3]
		y_modelo[4] = y_ini[4]
		y_modelo[5] = y_ini[5]
		y_modelo[6] = y_ini[6]
		
		Ca = y_modelo[1]
		Cb = y_modelo[2]
		Cc = y_modelo[3]
		Cd = y_modelo[4]
		T  = y_modelo[5]
		Tc = y_modelo[6]
		
		*/
		
	CONTINUOUS
	
		k10 = k10_estim*1e9
		k20 = k20_estim*1e10
		k30 = k30_estim*1e11
		-- Lectura de variables
		-- Variables de entrada
		q = linearInterp1D(tab_q, TIME)
		Fr = linearInterp1D(tab_Fr, TIME)
		T0 = linearInterp1D(tab_T0, TIME)
		Tc0 = linearInterp1D(tab_Tc0, TIME)
		-- Variables de salida
		y_real[1] = linearInterp1D(tab_Ca, TIME)
		y_real[2] = linearInterp1D(tab_Cb, TIME)
		y_real[3] = linearInterp1D(tab_Cc, TIME)
		y_real[4] = linearInterp1D(tab_Cd, TIME)
		y_real[5] = linearInterp1D(tab_T, TIME)
		y_real[6] = linearInterp1D(tab_Tc, TIME)		
		
		-- Balances de materia 
		V*Ca' = q*(Ca0 - Ca) + V*(-r1 - 2*r3)
		V*Cb' = -q*Cb + V*(r1 - r2)
		V*Cc' = -q*Cc + V*r2
		V*Cd' = -q*Cd + V*r3
		
		r1 = k10*exp(-Ea1/(R*(T+273.15)))*Ca
		r2 = k20*exp(-Ea2/(R*(T+273.15)))*Cb
		r3 = k30*exp(-Ea3/(R*(T+273.15)))*Ca**2
		
		-- Balances de energía
		-- Reactor
		rho*Cp*V*T' = q*rho*Cp*(T0 - T) - Q + Heat_rxn	
		-- Camisa
		rho*Cp*Vc*Tc' = Fr*rho*Cp*(Tc0 - Tc) + Q
		-- Transferencia de calor
		Q = UA*(T - Tc)
		--UA = alpha
		UA = alpha*(Fr**0.8)
		-- Calor generado por la reacción
		Heat_rxn = V*( - dHrxn1*r1 - dHrxn2*r2 - 2*dHrxn3*r3 )
		
		-- Optimizador	
		-- Subtotales del índice de coste
		EXPAND(i IN 1, nsal) J[i] = (abs(y_modelo[i] - y_real[i])/(media[i]*C) - log(1+abs(y_modelo[i] - y_real[i])/(media[i]*C)))
		-- Funcion de coste
      J_costo = C**2*SUM(i IN 1, nsal; pesos[i]*J[i])
		
		-- Restricciones no lineales (expresiones  c(x) <= 0)
		F_optim[1] = J_costo
		F_optim[2] = MATH.min(T,LiminfT) - LiminfT
	   F_optim[3] = MATH.max(T,LimsupT) - LimsupT
	   F_optim[4] = MATH.min(Tc,LiminfTc) - LiminfTc
	   F_optim[5] = MATH.max(Tc,LimsupTc) - LimsupTc
	   F_optim[6] = MATH.min(Ca,LiminfCa) - LiminfCa
	   F_optim[7] = MATH.max(Ca,LimsupCa) - LimsupCa	
		
		-- Asignación de valores del modelo
		y_modelo[1] = Ca
		y_modelo[2] = Cb
		y_modelo[3] = Cc
		y_modelo[4] = Cd
		y_modelo[5] = T
		y_modelo[6] = Tc


END COMPONENT