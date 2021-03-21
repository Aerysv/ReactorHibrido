/*-----------------------------------------------------------------------------------------
 LIBRARY: BLOQUECaLCULADOR
 FILE: Bloque_CalcuLectura
 AUTHOR: Daniel Montes
 COMPANY: Universidad de Valladolid
 DESCRIPTION: Errores de dHrxn corregidos.
 CREATION DATE: 21/03/2021
 
 Existen dos modos de operación:
  * Modo 0: A -> B
  * Mood 1: A -> B -> C
  			   2A -> D
-----------------------------------------------------------------------------------------*/
COMPONENT Bloque_CalcuLectura

	DATA

		REAL V = 11.5					"Volumen del reactor(L)"
		REAL R = 0.00831				"Constante de los gases ideales (kJ/mol/K)"
		REAL Ca0 = 5					"Concentración entrada A (mol/L)"
		-- Parámetros cinéticos
		REAL k10 = 6.4e9				"Factor preexponencial (1/min)"
		REAL k20 = 4.84e10			"Factor preexponencial (1/min)"
		REAL k30 = 8.86e11			"Factor preexponencial (1/min/mol)"
		REAL Ea1 = 52.1				"Energía de activacion kJ/mol"
		REAL Ea2 = 70.0				"Energía de activacion kJ/mol"
		REAL Ea3 = 65.0				"Energía de activacion kJ/mol"
		REAL dHrxn1 = -28.75			"Calor de reaccion (kJ/mol)"
		REAL dHrxn2 = -34.50			"Calor de reaccion (kJ/mol)"
		REAL dHrxn3 = -40.25			"Calor de reaccion (kJ/mol)"
		
		REAL q = 1
		REAL T = 20
		REAL Operacion = 1
	
	DECLS
	
	-- Calculos
		REAL Ca, Cb, Cc, Cd			"Concentración A (mol/L)"
		REAL r1, r2, r3				"Velocidad de la reaccion 1 (mol/L/min)"		
	-- Salidas del bloque Calculador	
		REAL x		"conversion (s/u)"	
		REAL P		"Potencia proporcionada por la bomba = Calor de reacción (W)"
				
	INIT
		
		Ca = 0
		Cb = 0
		Cc = 0
		Cd = 0

	CONTINUOUS
	-- Las reacciones son
	-- La variable Operacion "apaga" o "enciende" 2 y 3 para simular Van de Vusse.
	-- Calculo de la conversion
		x = (Ca0 - Ca)/Ca0
	-- Calculo de la velocidad de reacción
		r1 = k10*exp(-Ea1/(R*(T+273.15)))*Ca
		r2 = Operacion*k20*exp(-Ea2/(R*(T+273.15)))*Cb
		r3 = Operacion*k30*exp(-Ea3/(R*(T+273.15)))*Ca**2
	-- Balances de materia 
		V*Ca' = (q*(Ca0 - Ca) + V*(-r1 - 2*r3))/60
		V*Cb' = (-q*Cb + V*( r1 -   r2))/60
		V*Cc' = (-q*Cc + V*( r2       ))/60
		V*Cd' = (-q*Cd + V*(        r3))/60
	-- Potencia que proporciona la resistencia
		-- La potencia debe salir en vatios pero el término V*dHrxni*ri tiene unidades L*kJ/mol*mol/L/min = kJ/min
		-- Se apliCa *1000 para convertir a J y se divide por 60 para convertir a segundos. J/s  = W
		P = (1000/60)*V*( - dHrxn1*r1 - dHrxn2*r2 - 2*dHrxn3*r3 )
				
END COMPONENT