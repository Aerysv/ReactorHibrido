COMPONENT Reactor_RTO

	DATA
		-- Parametros del modelo
		REAL V = 11.5			"Volumen del reactor (L)"
		REAL Vc = 1				"Volumen de la camisa (L)"
		REAL R = 0.00831		"Constante de los gases (kJ/mol/K)"
		REAL rho = 1			"Densidad del agua (kg/L)"
		REAL Cp = 4.184		"Capacidad Calorifica del agua (kJ/kg/C)"	
      REAL Ca0 = 5			"Concentracion de entrada de A (mol/L)"
		REAL T0 = 15			"Temperatura de entrada de los reactivos (C)"
		REAL Tc0 = 15			"Temperatura de entrada de los reactivos (C)"

		-- Parametros estimados
		REAL k10 = 9.94755854e+10	"Factor pre-exponencial reacción 1 (1/min)"
		REAL k20 = 9.98553963e+11	"Factor pre-exponencial reacción 1 (1/min)"
		REAL k30 = 9.99263023e+12	"Factor pre-exponencial reacción 1 (1/min/mol)"
		REAL Ea1 = 58.9922611		"Energía de activación reacción 1 (kJ/mol)"
		REAL Ea2 = 77.6157046		"Energía de activación reacción 2 (kJ/mol)"
		REAL Ea3 = 71.1106314		"Energía de activación reacción 3 (kJ/mol)"
		REAL dHrxn1 = -21.2199341	"Entalpía de reacción 1 (kJ/mol)"
		REAL dHrxn2 = -2.68152145	"Entalpía de reacción 1 (kJ/mol)"
		REAL dHrxn3 = -66.5367189	"Entalpía de reacción 1 (kJ/mol)"
		REAL alpha = 1.59324467		"Parámetro de transferencia de calor"

		-- Optimizador
		REAL LiminfT = 5
    	REAL LimsupT = 60
	 	REAL LiminfCb = 0
    	REAL LimsupCb = 5
		REAL Liminfq = 0.3
		REAL Limsupq = 1.2
		REAL LiminfFr = 4
		REAL LimsupFr = 15
		-- Variables de decisión		
		REAL q = 1.2
		REAL Fr = 7
		-- Costos
		REAL pA = 0.5
		REAL pB = 5
		REAL pFr = 0.3
	
	DECLS
		-- Variables de decisión
		REAL T = 30			"Temperatura del reactor (C)"
		REAL Tc  = 20     "Temperatura refrigerante (C)"
		REAL Ca = 0.2		"Concentración de A (mol/L)"
		REAL Cb = 3			"Concentración de B (mol/L)"
		REAL Cc = 0.3		"Concentración de C (mol/L)"
		REAL Cd = 0.2		"Concentración de D (mol/L)"
		REAL r1				"Velocidad de la reacción 1 (1/min)"
		REAL r2				"Velocidad de la reacción 2 (1/min)"
		REAL r3				"Velocidad de la reacción 3 (mol/L/min)"
		REAL Heat_rxn		"Calor generado por la reacción (kJ/min)"
		REAL Q				"Calor transferido con la camisa (kJ/min)"
		REAL UA				"Coeficiente de transferencia de calor (kJ/min/C)"
		-- Optimizador
		REAL F_optim[3]		"Restricciones"
		REAL J_costo			"Funcion de coste"

	INIT

	CONTINUOUS

	-- Balances de materia 
		0 = q*(Ca0 - Ca) + V*(-r1 - 2*r3)
		0 = -q*Cb + V*(r1 - r2)
		0 = -q*Cc + V*r2
		0 = -q*Cd + V*r3
		
		r1 = k10*exp(-Ea1/(R*(T+273.15)))*Ca
		r2 = k20*exp(-Ea2/(R*(T+273.15)))*Cb
		r3 = k30*exp(-Ea3/(R*(T+273.15)))*Ca**2

	-- Balances de energia
		-- Reactor
		0 = q*(T0 - T) + (- Q + Heat_rxn)/(rho*Cp)
		-- Camisa
		0 = Fr*(Tc0 - Tc) + Q/(rho*Cp)
		-- Transferencia de calor
		Q = UA*(T - Tc)
		UA = alpha*Fr**0.8
		-- Calor generado por la reaccion
		Heat_rxn = V*( - dHrxn1*r1 - dHrxn2*r2 - 2*dHrxn3*r3 )

	-- Esta funcion de costo es solo para visualizacion
		J_costo = -( q*(pB*Cb - pA*Ca0) - pFr*Fr )

		-- Restricciones no lineales (expresiones  c(x) <= 0)
		F_optim[1] = J_costo
		F_optim[2] = T
	   F_optim[3] = Cb

END COMPONENT
