/*-----------------------------------------------------------------------------------------
 LIBRARY: PLANTAHIBRIDA
 FILE: CostoMHE
 CREATION DATE: 17/07/2020
-----------------------------------------------------------------------------------------*/
USE MATH

COMPONENT CostoMHE (INTEGER Nx = 4, INTEGER Ne = 4, INTEGER Nm = 4, INTEGER Nd = 2)

-- Nu = Control horizon (number of sampling periods)
-- MV = Number of manipulated variables of the MPC
-- PV = Number of process variables with Upper/lower constraints
-- Nx = Number of states of the model
-- Nm = Number of measured process variables
-- Nd = Number of measured disturbances
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
	REAL Cp = 4.184		"Capacidad calorifica del agua (kJ/kg/C)"	
   REAL Ca0 = 5			"Concentracion de entrada de A (mol/L)"
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
	-- Parametros del MHE
	REAL t_Sample = 0.5				"Tiempo de muestro (min)"
	REAL uq[Ne] = 1.75
	REAL uFr[Ne] = 8.51
	REAL Sig = 8      "Parametro de la Sigmoide" 

	REAL Liminfq = 0.3
	REAL Limsupq = 2
	REAL LiminfFr = 0.83
	REAL LimsupFr = 15
	REAL LiminfT = 5
	REAL LimsupT = 60
	REAL LiminfCb = 0
	REAL LimsupCb = 5
	
	-- Tablas de datos
	REAL tt[5] = {0, 0.5, 1, 1.5, 2}    "Tiempos pasados base de table(min)"
	REAL ty[1] = 0 	
	REAL tz[1] = 0
	-- Medidas de las acciones de control
	REAL tq[Ne+1] = 1.2		"Valores pasados de medidas de q para tabla"
	REAL tFr[Ne+1] = 7		"Valores pasados de medidas de Fr para tabla"
	-- Medidas de las perturbaciones medidas
	REAL tT0[Ne+1] = 10		"Valores pasados de medidas de T0 para tabla"
	REAL tTc0[Ne+1] = 10		"Valores pasados de medidas de Tc0 para tabla"
	-- Estados medidos
	REAL tCa[Ne+1]= 0.55		"Valores pasados de medidas de Ca para tabla"
	REAL tCb[Ne+1]= 0.44		"Valores pasados de medidas de Cb para tabla"
	REAL tT[Ne+1]= 45			"Valores pasados de medidas de T para tabla"
	REAL tTc[Ne+1]= 15.59	"Valores pasados de medidas de Tc para tabla"
				
	REAL beta_xv = 0.1		"Importancia de las perturbaciones en el estimador"
	REAL beta_xN = 0.1		"Importancia de las variables medidas en el estimador"
	
DECLS

	REAL q               "Caudal reactivos (L/min)"
	REAL Fr              "Caudal refrigerante (L/min)"
	
	REAL T0					"Temperatura de entrada de reactivos (C)"
	REAL Tc0					"Temperatura del reactor (C)"
	
    -- Salidas medidas del proceso
	REAL T                 "Temperatura reactor (C)"
	REAL Tc                "Temperatura reFrigerante (C)"

	-- Otras variables del modelo
	REAL Q					"Calor transferido (kJ/mol)"
	REAL Heat_rxn			"Calor generado (kJ/min)"
	REAL Ca, Cb
	REAL r1, r2, r3
	REAL UA
	
-- Variables del problema de optimizacion
	REAL J_costo			"Costo MHE"
	REAL F_optim[5]		"Restricciones"	
	REAL J_costo_m 
	REAL J_costo_N 
	REAL J_costo_v 	

-- estado en t-N	que fue estimado en t-1

	DISCR REAL Ca_Ne = 0.55
	DISCR REAL Cb_Ne = 0.44
	DISCR REAL T_Ne = 45
	DISCR REAL Tc_Ne = 15.59

-- estado a estimar en t-N

	DISCR REAL Ca_N = 0.55
	DISCR REAL Cb_N = 0.44
	DISCR REAL T_N = 45
	DISCR REAL Tc_N = 15.59
	
-- valores medidos de la salida (estado) en t

	REAL Ca_m = 0.55
	REAL Cb_m = 0.44
	REAL T_m = 45
	REAL Tc_m = 15.59
	
	--  perturbaciones del estado a estimar 

	DISCR REAL v[Nx] = 0		    	"perturbaciones del estado a estimar"   
	DISCR REAL v_ant[Nx] = 0    	"perturbaciones del estado a estimar"  
	REAL  error[2]						"Error de prediccion de las salidas"
	
OBJECTS
	-- Acciones de control
	TABLE tMedq
	TABLE tMedFr
	-- Perturbaciones medidas
	TABLE tMedT0
	TABLE tMedTc0
	-- Estados medidos
	TABLE tMeda
	TABLE tMedb
	TABLE tMedT
	TABLE tMedTc	
	
INIT

		J_costo = 0
		J_costo_m = 0
		J_costo_N = 0
		J_costo_v = 0
		
		FOR (i IN 1,5)
			F_optim[i] = 0
		END FOR
		
    	Ca = Ca_N    
    	Cb = Cb_N       
    	T = T_N    
    	Tc = Tc_N 

--   4 = Ne +1		
		-- Medidas de las acciones de control
		tMedq.fillData("tableq",5,tt,0,ty,0,tz,5,tq)
		tMedFr.fillData("tableFr",5,tt,0,ty,0,tz,5,tFr)
		-- Perturbaciones medidas
		tMedT0.fillData("tableT0",5,tt,0,ty,0,tz,5,tT0)
		tMedTc0.fillData("tableTc0",5,tt,0,ty,0,tz,5,tTc0)		
		-- Estados medidos
		tMeda.fillData("tableCa",5,tt,0,ty,0,tz,5,tCa)
	 	tMedb.fillData("tableCb",5,tt,0,ty,0,tz,5,tCb)
	 	tMedT.fillData("tableT",5,tt,0,ty,0,tz,5,tT)
	 	tMedTc.fillData("tableTc",5,tt,0,ty,0,tz,5,tTc)
		
		
CONTINUOUS
		-- Medidas de las acciones de control
		q = tMedq.interp1D(CONSTANT,CONSTANT,TIME)
		Fr = tMedFr.interp1D(CONSTANT,CONSTANT,TIME)
		-- Perturbaciones medidas
		T0 = tMedT0.interp1D(SPLINE,SPLINE,TIME)
		Tc0 = tMedTc0.interp1D(SPLINE,SPLINE,TIME)	 
		-- Estados medidos
   	Ca_m = tMeda.interp1D(SPLINE,SPLINE,TIME)	
    	Cb_m = tMedb.interp1D(SPLINE,SPLINE,TIME)		
    	T_m =  tMedT.interp1D(SPLINE,SPLINE,TIME)	
    	Tc_m = tMedTc.interp1D(SPLINE,SPLINE,TIME)	
		
		-- Balances de materia 
		V*Ca' = q*(Ca0 - Ca) + V*(-r1 - 2*r3) + v[1]
		V*Cb' = -q*Cb + V*( r1 -   r2) + v[2]
		
		r1 = k10*exp(-Ea1/(R*(T+273.15)))*Ca
		r2 = k20*exp(-Ea2/(R*(T+273.15)))*Cb
		r3 = k30*exp(-Ea3/(R*(T+273.15)))*Ca**2
		
		-- Balances de energia
		-- Reactor
		V*T' = q*(T0 - T) + (- Q + Heat_rxn)/(rho*Cp) + v[3]
		-- Camisa
		Vc*Tc' = Fr*(Tc0 - Tc) + Q/(rho*Cp) + v[4]
		-- Transferencia de calor
		Q = UA*(T - Tc)
		UA = alpha*Fr**0.8
		-- Calor generado por la reaccion
		Heat_rxn = V*( - dHrxn1*r1 - dHrxn2*r2 - 2*dHrxn3*r3 )
		
-- Costo MHE

		J_costo_m =  ( (Ca/Ca_m -1)**2 + (Cb/Cb_m -1)**2 + \
						 (T/T_m -1)**2 + (Tc/Tc -1)**2 )  
		--J_costo_m = ((Ca_m - Ca)/(5))**2 + ((Cb_m - Cb)/(5))**2 + \
		--				((T_m - T)/(40))**2 + ((Tc_m - Tc)/(10))**2
		J_costo_N =  (Ca_N/Ca_Ne -1)**2 + (Cb_N/Cb_Ne -1)**2 + \
			 	    (T_N/T_Ne -1)**2 + (Tc_N/Tc_Ne -1)**2
		--J_costo_N = 0
		--J_costo_v =  SUM(i IN 1,Nx; (v[i] - v_ant[i])**2)
		J_costo_v =  SUM(i IN 1, Nx; (v[i])**2)
		
		J_costo =   100*J_costo_m + beta_xN* J_costo_N + beta_xv*J_costo_v  
	
		error[1] =  T_m - T
		error[2] =  Cb_m - Cb
	
		-- Restricciones no lineales (expresiones  c(x) <= 0)
		F_optim[1] = J_costo
		F_optim[2] = MATH.min(T, LiminfT) - LiminfT
	   F_optim[3] = MATH.max(T, LimsupT) - LimsupT
	   F_optim[4] = MATH.min(Cb, LiminfCb) - LiminfCb
	   F_optim[5] = MATH.max(Cb, LimsupCb) - LimsupCb	
	
END COMPONENT