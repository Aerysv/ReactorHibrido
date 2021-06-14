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
	REAL Pred_h = 60		"Horizonte de predicción (min)"
		
	-- Precios del controlador
	REAL T_sp = 35
	REAL Cb_sp = 2.6
	
	-- Variables manipuladas
	REAL beta[2] =  {0.5, 0.5}		"Penalizacion de cambios"
	REAL gamma[2] = {10, 10}				"Importancia relativa de los setpoint"
	
	-- Parámetros estimador
	REAL margen_v = 15		"Rango de variación de las perturbaciones"
	REAL margen_x = 20		"Rango de variación del estado inicial del MHE (%)"	
	REAL beta_xv = 0.01		"Peso de las perturbaciones en el costo del MHE"
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
	DISCR REAL T = 35                 "Temperatura reactor (C)"
	DISCR REAL Tc = 25               "Temperatura reFrigerante (C)"

	-- Otras variables del modelo
	DISCR REAL Ca = 0.05
	DISCR REAL Cb = 2.5
	
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
	DISCR REAL state[4] = {0.05, 2.5, 35, 25}				--  estimated (state) variables (+ algebraic)
	DISCR REAL u_new[MV*Nu]      	--  Vector of decision variables computed by the MPC + slacks
	DISCR REAL v_new[Nx] = 0   			--  valores estimados de las perturbaciones
	DISCR REAL acc[MV] 				-- Valores actuales de las acciones de control
	DISCR REAL per[Nd]				-- Valores actuales de las perturbaciones medidas
	DISCR REAL med[Nm]				-- current measurements
	DISCR REAL acc_ant[MV*(Ne+1)]		--  Acciones de control medidas pasadas
	DISCR REAL per_ant[Nd*(Ne+1)]		--  Perturbaciones medidas pasadas
	DISCR REAL med_ant[Nm*(Ne+1)]		--  past measured variables
	DISCR REAL aux[4]							-- valores auxiliares

	DISCR REAL u_ant[MV*Ne]				-- valores pasados del control aplicado
	DISCR REAL x_ant[Nx,Ne]				--  valores de los estados del modelo en el periodo anterior
	DISCR REAL x_Ne[4] = {0.05, 2.5, 35, 25}					-- Estimated state at t- Ne
	
-- Upper and lower bounds of the manipulated and controlled variables	
	DISCR REAL lim_manip_up[MV], lim_manip_low[MV], lim_con_up[PV], lim_con_low[PV]
	DISCR REAL uqant = 0.9				"Accion de control anterior"
	DISCR REAL uFrant = 0				"Accion de control anterior"
	DISCR REAL config[6]			"Configuracion del controlador"
	DISCR REAL Pred_hh = 60		"Horizonte de predicción (min)"
	DISCR REAL error[2] = 0
	
-- Acciones de control
	DISCR REAL uq[Nu] = 0.9	"Acciones de control de reactivos (mol/L)"
	DISCR REAL uFr[Nu] = 9	"Acciones de control de refrigerante (mol/L)"
	
OBJECTS 

	optimizCON modeloCON
	estimMHE   estadoMHE
		
END COMPONENT