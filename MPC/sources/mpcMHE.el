/*-----------------------------------------------------------------------------------------
 LIBRARY: PLANTAHIBRIDA
 FILE: mpcMHE
 CREATION DATE: 09/07/2020
-----------------------------------------------------------------------------------------*/
USE MATH

COMPONENT mpcMHE (INTEGER Nu = 3, INTEGER MV = 2, INTEGER PV = 2, INTEGER Nx = 4, INTEGER Nm = 4, INTEGER Nd = 2, INTEGER Ne = 4)
														

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
-- B ==> C
-- 2A ==> D
-- Cb es el producto deseado
-- Se miden todas las variables
-- Modelo perfecto


DATA

	-- Parametros del modelo
	REAL V = 11.5			"Volumen del reactor (L)"
	REAL Vc = 1				"Volumen de la camisa (L)"
	REAL R = 0.00831		"Constante de los gases (kJ/mol/K)"
	REAL rho = 1			"Densidad del agua (kg/L)"
	REAL Cp = 4.184		"Capacidad Calorifica del agua (kJ/kg/C)"	
   REAL Ca0 = 5			"Concentracion de entrada de A (mol/L)"
	REAL T0 = 10			"Temperatura de entrada de los reactivos (C)"
	REAL Tc0 = 10			"Temperatura de entrada de los reactivos (C)"
	
	-- Parametros usados
	
	REAL k10 = 6.4e9 "Factor preexponencial (1/min)"
	REAL k20 = 4.84e10 "Factor preexponencial (1/min)"
	REAL k30 = 8.86e11 "Factor preexponencial (1/min/mol)"
	REAL Ea1 = 52.1 "Energía de activacion kJ/mol"
	REAL Ea2 = 70.0 "Energía de activacion kJ/mol"
	REAL Ea3 = 65.0 "Energía de activacion kJ/mol"
	REAL dHrxn1 = -28.75 "Calor de reaccion (kJ/mol)"
	REAL dHrxn2 = -34.50 "Calor de reaccion (kJ/mol)"
	REAL dHrxn3 = -40.25 "Calor de reaccion (kJ/mol)"
	REAL alpha = 1.8

	-- Parametros del controlador
	REAL t_Sample = 0.5				"Tiempo de muestro (min)"
	REAL Pred_h = 30
	
	REAL Liminfq = 0.3
	REAL Limsupq = 2
	REAL LiminfFr = 0.83
	REAL LimsupFr = 15
	REAL LiminfT = 5
	REAL LimsupT = 60
	REAL LiminfCb = 0
	REAL LimsupCb = 5
	
	-- Consignas
	REAL T_sp = 30
	REAL Cb_sp = 3.3
	REAL gamma[2] = {10, 10}				"Importancia relativa de los setpoint"
	-- Variables manipuladas
	REAL beta[2] =  {0.5, 0.5}				"Penalizacion de cambios"
	
	-- Parámetros estimador
	REAL margen_v = 0.5       	   "Rango de variación de las perturbaciones"
	REAL margen_x = 20				"Rango de variación del estado inicial del MHE (%)"	
	REAL beta_xv = 0.1			   "Peso de las perturbaciones en el costo del MHE"
	REAL beta_xN = 1					"Peso del estado en t-N en el costo de MHE"
	REAL v_ini = 0						"Inicializacion vector de perturbaciones"
	
DECLS

	-- Caudales
	REAL q                "Caudal reactivos (L/min)"
	REAL Fr               "Caudal reFrigerante (L/min)"
	
    -- Salidas medidas del proceso
	REAL T                 "Temperatura reactor (C)"
	REAL Tc                "Temperatura reFrigerante (C)"

	-- Otras variables del modelo
	REAL Q					"Calor transferido (kJ/mol)"
	REAL Heat_rxn			"Calor generado (kJ/min)"
	REAL Ca, Cb				"Concentraciones de las especies (mol/L)"
	REAL r1, r2, r3		"Velocidades de reacción (1/min)"
	REAL UA					"Coeficiente de transmisión de calor (kJ/C)"
	
-- Variables del problema de optimizacion

	BOOLEAN sample = TRUE
	BOOLEAN ControlFlag = FALSE

	DISCR REAL state[Nx]					--  estimated (state) variables (+ algebraic)
	DISCR REAL u_new[MV*Nu]      		--  Vector of decision variables computed by the MPC + slacks
	DISCR REAL v_new[Nx]   				--  valores estimados de las perturbaciones
	REAL acc[MV]							-- Valores actuales de las acciones de control
	REAL per[Nd]							-- Valores actuales de las perturbaciones medidas
	REAL med[Nm]							-- current measurements
	DISCR REAL acc_ant[MV*(Ne+1)]		--  Acciones de control medidas pasadas
	DISCR REAL per_ant[Nd*(Ne+1)]		--  Perturbaciones medidas pasadas
	DISCR REAL med_ant[Nm*(Ne+1)]		--  past measured variables
	DISCR REAL aux[4]							-- valores auxiliares

	DISCR REAL u_ant[MV*Ne]				-- valores pasados del control aplicado
	DISCR REAL x_ant[Nx,Ne]				--  valores de los estados del modelo en el periodo anterior
	DISCR REAL x_Ne[Nx]					-- Estimated state at t- Ne
	
-- Upper and lower bounds of the manipulated and controlled variables	
	DISCR REAL lim_manip_up[MV], lim_manip_low[MV], lim_con_up[PV], lim_con_low[PV]
	DISCR REAL uqant				"Accion de control anterior"
	DISCR REAL uFrant				"Accion de control anterior"
	DISCR REAL config[6]			"Configuracion del controlador"
	DISCR REAL Pred_hh = 30
	DISCR REAL error[2]
	
-- Acciones de control
	DISCR REAL uq[Nu] = 1.75
	DISCR REAL uFr[Nu] = 8.51
	
OBJECTS 

	optimizCON modeloCON
	estimMHE   estadoMHE
		
INIT

-- Condiciones iniciales de las variables de estado
      Ca = 0.08
      Cb = 2.3
      T = 25
      Tc = 13

-- constraints MVs
  		lim_manip_low[1] = Liminfq        -- limites  q
		lim_manip_up[1] = Limsupq

		lim_manip_low[2] = LiminfFr		 -- limites  Fr		
  	 	lim_manip_up[2] = LimsupFr

-- Constraints process variables
		lim_con_up[1]  = LimsupT          -- limites temperatura reactor
		lim_con_low[1] = LiminfT
	
		lim_con_up[2] = LimsupCb           -- limites Cb
		lim_con_low[2] = LiminfCb
	
-- Initial values MVs

	FOR (i IN 1,Nu)
		u_new[MV*(i-1)+1] = uq[i]
		u_new[MV*(i-1)+2] = uFr[i]
	END FOR
	
	FOR (j IN 1,Ne)
			u_ant[MV*(j-1)+1] = uq[1]
			u_ant[MV*(j-1)+2] = uFr[1]			
	END FOR	

		state[1] = Ca
		state[2] = Cb										
		state[3] = T
		state[4] = Tc
		
		acc[1] = q
		acc[2] = Fr
		
		per[1] = T0
		per[2] = Tc0
		
		med[1] = Ca
		med[2] = Cb									
		med[3] = T
		med[4] = Tc	
	
		config[1] = T_sp
		config[2] = Cb_sp
		config[3] = gamma[1]
		config[4] = gamma[2]
		config[5] = beta[1]
		config[6] = beta[2]
		
		q = uq[1]
		Fr = uFr[1]

		FOR (j IN 1,Ne+1)		
			acc_ant[MV*(j-1)+1] = q
			acc_ant[MV*(j-1)+2] = Fr
			
			per_ant[Nd*(j-1)+1] = T0
			per_ant[Nd*(j-1)+2] = Tc0

			med_ant[Nm*(j-1)+1] = Ca		
			med_ant[Nm*(j-1)+2] = Cb		
			med_ant[Nm*(j-1)+3] = T		
			med_ant[Nm*(j-1)+4] = Tc					
		END FOR		

		FOR (j IN 1,Ne)		
			FOR (i IN 1,Nx)
				x_ant[i,j] = state[i]
			END FOR		
		END FOR
		
		FOR (i IN 1,Nx)
			x_Ne[i] = x_ant[i,1]
			v_new[i] = v_ini
		END FOR		
		
		FOR (i IN 1,2)
			error[i] = 0
		END FOR		
		
-- Initialization of the sensitivities	

   		modeloCON.iniciar()		
			estadoMHE.iniciar()
		
DISCRETE

	-- Llamada al controlador predictivo

	WHEN (sample == TRUE) THEN

		uqant = q
		uFrant= Fr

		aux[1] = uqant
		aux[2] = uFrant
		aux[3] = T0
		aux[4] = Tc0

		-- Actualizar el vector de acciones de control
		FOR(i IN 1, MV)
			FOR (j IN 1,Ne)
				acc_ant[MV*(j-1)+i] = acc_ant[MV*j+i]
			END FOR		
			acc_ant[MV*Ne+i] = acc[i]
		END FOR			
		
		-- Actualizar el vector de perturbaciones medidas
		FOR(i IN 1, Nd)
			FOR (j IN 1,Ne)
				per_ant[Nd*(j-1)+i] = per_ant[Nd*j+i]
			END FOR		
			per_ant[Nd*Ne+i] = per[i]
		END FOR		
		
		-- actualizar vector de medidas
		FOR (i IN 1,Nm)
			FOR (j IN 1,Ne)
				med_ant[Nm*(j-1)+i] = med_ant[Nm*j+i]
			END FOR
			med_ant[Nm*Ne+i] = med[i]
		END FOR	
		
-- actualizar estados anteriores			
/*
		state[1] = Ca
		state[2] = Cb
		state[3] = T
		state[4] = Tc
		*/
		FOR (i IN 1,Nx)
			FOR (j IN 1,Ne-1)
				x_ant[i,j] = x_ant[i,j+1]
			END FOR
			x_ant[i,Ne] = state[i]
		END FOR	
		
		FOR (i IN 1,Nx)
				x_Ne[i] = x_ant[i,1]
		END FOR	
		
-- 	actualizar controles anteriores aplicados	

		FOR (i IN 1,MV)
			FOR (j IN 1,Ne-1)
				u_ant[MV*(j-1)+i] = u_ant[MV*j+i]
			END FOR
				u_ant[MV*(Ne-1)+i] = u_new[i]
		END FOR
				
		sample = FALSE
		sample = TRUE AFTER t_Sample
		ControlFlag = TRUE AFTER (Ne*t_Sample)

	END WHEN
	
	WHEN (ControlFlag == TRUE) THEN
	
		config[1] = T_sp
		config[2] = Cb_sp
		config[3] = gamma[1]
		config[4] = gamma[2]
		config[5] = beta[1]
		config[6] = beta[2]	
		
		Pred_hh = Pred_h
		
		-- constraints MVs
   	lim_manip_low[1] = Liminfq        -- limites  q
		lim_manip_up[1] = Limsupq

		lim_manip_low[2] = LiminfFr		 -- limites  Fr		
   	lim_manip_up[2] = LimsupFr

-- Constraints process variables
		lim_con_up[1]  = LimsupT          -- limites temperatura reactor
		lim_con_low[1] = LiminfT
	
		lim_con_up[2] = LimsupCb           -- limites Cb
		lim_con_low[2] = LiminfCb
		
		-- Estimar y recoger los estados y perturbaciones
		estadoMHE.mhe (t_Sample, x_Ne, acc_ant, per_ant, med_ant,  u_ant,  margen_v,  margen_x,  beta_xv, beta_xN,  state, v_new, error )
				
		-- Calcular y recoger la nueva señal de control
		modeloCON.control (t_Sample, Pred_hh, state, aux, v_new, error, lim_manip_low, lim_manip_up, lim_con_low, lim_con_up, config, u_new )

		FOR (i IN 1,Nu)
			uq[i] = u_new[MV*(i-1)+1]
			uFr[i] = u_new[MV*(i-1)+2]
		END FOR	
		
		ControlFlag = FALSE
	END WHEN
	
CONTINUOUS

		-- Actualizacion acciones de control
		q = uq[1]
		Fr = uFr[1]
		
		-- Balances de materia 
		V*Ca' = q*(Ca0 - Ca) + V*(-r1 - 2*r3)
		V*Cb' = -q*Cb + V*( r1 -   r2)
		
		r1 = k10*exp(-Ea1/(R*(T+273.15)))*Ca
		r2 = k20*exp(-Ea2/(R*(T+273.15)))*Cb
		r3 = k30*exp(-Ea3/(R*(T+273.15)))*Ca**2
		
		-- Balances de energia
		-- Reactor
		rho*Cp*V*T' = q*rho*Cp*(T0 - T) - Q + Heat_rxn	
		-- Camisa
		rho*Cp*Vc*Tc' = Fr*rho*Cp*(Tc0 - Tc) + Q
		-- Transferencia de calor
		Q = UA*(T - Tc)
		UA = alpha*Fr**0.8
		-- Calor generado por la reaccion
		Heat_rxn = V*( - dHrxn1*r1 - dHrxn2*r2 - 2*dHrxn3*r3 )
		
		acc[1] = q
		acc[2] = Fr
		
		per[1] = T0
		per[2] = Tc0
		
		med[1] = Ca 
		med[2] = Cb									
		med[3] = T
		med[4] = Tc	

	
END COMPONENT