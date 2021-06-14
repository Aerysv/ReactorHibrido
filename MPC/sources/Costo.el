/*-----------------------------------------------------------------------------------------
 LIBRARY: PLANTAHIBRIDA
 FILE: Costo
 CREATION DATE: 07/07/2020
-----------------------------------------------------------------------------------------*/
USE MATH

COMPONENT Costo (INTEGER Nu = 3, INTEGER MV = 2, INTEGER Nx = 4)

-- Nu = Control horizon (number of sampling periods)
-- MV = Number of manipulated variables of the MPC
-- PV = Number of process variables with Upper/lower constraints
-- Nx = Number of states of the model
-- Nm = Number of measured process variables 
-- Ns = Number of constrints with slack variables for feasibility 
-- Ne = Number of (past and current) samples used in MHE

-- Proceso con tres reacciones
-- A  ==> B
--  B ==> C
-- 2A ==> D
-- B es el producto deseado
-- Se miden las concentraciones y temperaturas
-- Modelo perfecto

DATA

-- Parametros del modelo
	REAL V = 11.5			"Volumen del reactor (L)"
	REAL Vc = 1				"Volumen de la camisa (L)"
	REAL R = 0.00831		"Constante de los gases (kJ/mol/K)"
	REAL rho = 1			"Densidad del agua (kg/L)"
	REAL Cp = 4.184		"Capacidad Calorifica del agua (kJ/kg/C)"	
   REAL Ca0 = 5			"Concentracion de entrada de A (mol/L)"
	REAL T0 = 23			"Temperatura de entrada de los reactivos (C)"
	REAL Tc0 = 23			"Temperatura de entrada de los reactivos (C)"
	/*
	-- Parametros estimados
	REAL k10 = 2.77278263e+10
	REAL k20 = 2.95155307e+11
	REAL k30 = 2.96918831e+12
	REAL Ea1 = 58.7887775
	REAL Ea2 = 75.335899
	REAL Ea3 = 73.8049428
	REAL dHrxn1 = -21.096637
	REAL dHrxn2 = -26.1608589
	REAL dHrxn3 = -47.3443562
	REAL alpha = 1.54979214
		*/
	REAL Ea1 = 58.9922611
	REAL Ea2 = 77.6157046
	REAL Ea3 = 71.1106314
	REAL alpha = 1.59324467
	REAL dHrxn1 = -21.2199341
	REAL dHrxn2 = -2.68152145
	REAL dHrxn3 = -66.5367189
	REAL k10 = 9.94755854e+10
	REAL k20 = 9.98553963e+11
	REAL k30 = 9.99263023e+12
	-- Parametros del controlador
	REAL T_sp = 35
	REAL Cb_sp = 2.6
	REAL gamma[2] = {1, 1}		"Importancia relativa de los setpoint"
	REAL beta[2] =  {0.1, 0.1}	"Penalizacion de cambios"
	REAL t_Sample = 0.5			"Tiempo de muestro (min)"
	REAL Liminfq = 0.3			"Límite inferior q (L/min)"
	REAL Limsupq = 1.2				"Límite superior q (L/min)"
	REAL LiminfFr = 5			"Límite inferior Fr (L/min)"
	REAL LimsupFr = 15			"Límite superior Fr (L/min)"
	REAL LiminfT = 5				"Límite inferior T (C)"
	REAL LimsupT = 60				"Límite superior T (C)"
	REAL LiminfCb = 0				"Límite inferior Cb (mol/L)"
	REAL LimsupCb = 5				"Límite superior Cb (mol/L)"
	
	REAL Sig = 100	"Parametro de la Sigmoide" 
	
DECLS

	REAL q			"Caudal reactivos (L/min)"
	REAL Fr 			"Caudal refrigerante (L/min)"
	
    -- Salidas medidas del proceso
	REAL T			"Temperatura reactor (C)"
	REAL Tc			"Temperatura reFrigerante (C)"

	-- Otras variables del modelo
	REAL Q			"Calor transferido (kJ/mol)"
	REAL Heat_rxn	"Calor generado (kJ/min)"
	REAL Ca, Cb		"Concentraciones (mol/L)"
	REAL r1, r2, r3		"Velocidades de reacción (1/min)"
	REAL UA			"Coeficiente de transmición de calor (kJ/C)"
	
-- Variables del problema de optimizacion
	REAL J_costo			"Función de costo"
	REAL F_optim[5]		"Restricciones"
	REAL J_setpoints		"Componente funcion de coste"
	REAL J_cambios			"Componente funcion de coste"
-- Variables del estimador de estados
	DISCR REAL v[Nx] = 0	"Perturbaciones estimadas"
	-- Estados iniciales estimador
	DISCR REAL Ca_inic = 0.05	"Concentración inicial de A (mol/L)"
	DISCR REAL Cb_inic = 2.6	"Concentración inicial de B (mol/L)"
	DISCR REAL T_inic = 35		"Temperatura inicial del reactor (C)"
	DISCR REAL Tc_inic = 28	"Temperatura inicial del refrigerante (C)"
	DISCR REAL error[2] = 0			"Error final del MHE"
-- Acciones de control
	DISCR REAL uq[Nu] = 0.9		"Caudal de reactivos (L/min)"
	DISCR REAL uFr[Nu] = 9		"Caudal de refrigerante (L/min)"
	DISCR REAL uqant, uFrant	"Acciones de control pasadas (L/min)"
	
INIT

	Ca = Ca_inic
	Cb = Cb_inic
	T = T_inic
	Tc = Tc_inic

CONTINUOUS
		-- Discretización de las acciones de control
		q = uq[1] + (uq[2]-uq[1])*sigmoide(TIME,t_Sample,Sig) + (uq[3]-uq[2])*sigmoide(TIME,2*t_Sample,Sig)
		Fr = uFr[1] + (uFr[2]-uFr[1])*sigmoide(TIME,t_Sample,Sig) + (uFr[3]-uFr[2])*sigmoide(TIME,2*t_Sample,Sig)
		
		V*Ca' = q*(Ca0 - Ca) + V*(-r1 - 2*r3) + v[1]
		V*Cb' = -q*Cb + V*(r1 - r2) + v[2]
		
		r1 = k10*exp(-Ea1/(R*(T+273.15)))*Ca
		r2 = k20*exp(-Ea2/(R*(T+273.15)))*Cb
		r3 = k30*exp(-Ea3/(R*(T+273.15)))*Ca**2
		
		-- Balances de energía
		-- Reactor
		rho*Cp*V*T' = q*rho*Cp*(T0 - T) - Q + Heat_rxn	+ v[3]
		-- Camisa
		rho*Cp*Vc*Tc' = Fr*rho*Cp*(Tc0 - Tc) + Q + v[4]
		-- Transferencia de calor
		Q = UA*(T - Tc)
		--UA = alpha
		UA = alpha*(Fr**0.8)
		-- Calor generado por la reacción
		Heat_rxn = V*( - dHrxn1*r1 - dHrxn2*r2 - 2*dHrxn3*r3 )
	
		J_setpoints = gamma[1]*((T + error[1] - T_sp)/(LimsupT-LiminfT))**2 + gamma[2]*((Cb + error[2] - Cb_sp)/(LimsupCb - LiminfCb))**2
		 
		J_cambios = beta[1]*((uq[1] - uqant)**2 + SUM(i IN 2,3; (uq[i]-uq[i-1])**2))/((Limsupq-Liminfq)**2) + \
					 + beta[2]*((uFr[1] - uFrant)**2 + SUM(i IN 2,3; (uFr[i]-uFr[i-1])**2))/((LimsupFr-LiminfFr)**2)
					 
		J_costo = J_setpoints + J_cambios
		
		-- Restricciones no lineales (expresiones  c(x) <= 0)
		F_optim[1] = J_costo
		F_optim[2] = MATH.min(T,LiminfT) - LiminfT
	   F_optim[3] = MATH.max(T,LimsupT) - LimsupT
	   F_optim[4] = MATH.min(Cb,LiminfCb) - LiminfCb
	   F_optim[5] = MATH.max(Cb,LimsupCb) - LimsupCb	
	
END COMPONENT