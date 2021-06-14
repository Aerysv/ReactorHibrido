USE MATH

COMPONENT mpcMHE_OPC (INTEGER Nu = 3, INTEGER MV = 2, INTEGER PV = 2, INTEGER Nx = 4, INTEGER Nm = 4, INTEGER Nd = 2, INTEGER Ne = 4)
														
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

DATA
	
	-- Parametros del controlador
	REAL t_Sample = 0.5	"Tiempo de muestro (min)"
	REAL Pred_h = 30		"Horizonte de predicción (min)"
		
	-- Precios del controlador
	REAL T_sp = 30
	REAL Cb_sp = 3.3
	
	-- Variables manipuladas
	REAL beta[2] =  {2, 2}		"Penalizacion de cambios"
	REAL gamma[2] = {10, 10}				"Importancia relativa de los setpoint"
	
	-- Parámetros estimador
	REAL margen_v = 0.5		"Rango de variación de las perturbaciones"
	REAL margen_x = 20		"Rango de variación del estado inicial del MHE (%)"	
	REAL beta_xv = 0.1		"Peso de las perturbaciones en el costo del MHE"
	REAL beta_xN = 1			"Peso del estado en t-N en el costo de MHE"
	REAL v_ini = 0				"Inicializacion vector de perturbaciones"
	
DECLS
	-- Perturbaciones
	DISCR REAL T0					"Temperatura de entrada de los reactivos (C)"
	DISCR REAL Tc0					"Temperatura de entrada de los refrigerante (C)"
	DISCR REAL Ca0					"Concentración de entrada de A (mol/L)"

	-- Caudales
	DISCR REAL q                "Caudal reactivos (L/min)"
	DISCR REAL Fr               "Caudal reFrigerante (L/min)"
	
    -- Salidas medidas del proceso
	DISCR REAL T                 "Temperatura reactor (C)"
	DISCR REAL Tc                "Temperatura reFrigerante (C)"

	-- Otras variables del modelo
	DISCR REAL Ca, Cb				"Concentraciones de las especies (mol/L)"
	
-- Variables del problema de optimizacion

	BOOLEAN sample = TRUE
	BOOLEAN ControlFlag = FALSE
	
	DISCR REAL Liminfq = 0.3		"Límite inferior de caudal de reactivos (L/min)"
	DISCR REAL Limsupq = 1.2		"Límite superior de caudal de reactivos (L/min)"
	DISCR REAL LiminfFr = 5		"Límite inferior de caudal de refrigerante (L/min)"
	DISCR REAL LimsupFr = 15		"Límite superior de caudal de refrigerante (L/min)"
	DISCR REAL LiminfT = 5			"Límite inferior de temperatura del reactor (L/min)"
	DISCR REAL LimsupT = 60			"Límite superior de temperatura del reactor (L/min)"
	DISCR REAL LiminfCb = 0			"Límite inferior de concentración de B (L/min)"
	DISCR REAL LimsupCb = 5			"Límite superior de concentración de B (L/min)"
	DISCR REAL state[Nx]				--  estimated (state) variables (+ algebraic)
	DISCR REAL u_new[MV*Nu]      	--  Vector of decision variables computed by the MPC + slacks
	DISCR REAL v_new[Nx]   			--  valores estimados de las perturbaciones
	DISCR REAL acc[MV]				-- Valores actuales de las acciones de control
	DISCR REAL per[Nd]				-- Valores actuales de las perturbaciones medidas
	DISCR REAL med[Nm]				-- current measurements
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
	DISCR REAL Pred_hh = 30		"Horizonte de predicción (min)"
	DISCR REAL error[2]
	
-- Acciones de control
	DISCR REAL uq[Nu] = 0.9	"Acciones de control de reactivos (mol/L)"
	DISCR REAL uFr[Nu] = 9	"Acciones de control de refrigerante (mol/L)"
	
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
	
END COMPONENT