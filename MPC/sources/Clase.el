/*-----------------------------------------------------------------------------------------
 LIBRARY: PLANTAHIBRIDA
 FILE: Clase
 CREATION DATE: 07/07/2020
-----------------------------------------------------------------------------------------*/
USE MATH
USE OPTIM_METHODS
		
CLASS optimizCON IS_A Costo_par1

	DECLS

		REAL Pred_h 
		CONST INTEGER numC = 4			-- número de restricciones del problema de optimización
		CONST INTEGER numU = 6			-- número de variables de decisión  MV*Nu + slacks
		INTEGER numX = 6			-- número de variables de estado del modelo	

	OBJECTS		
		
		EIDAS calculoSens

   METHODS

-- llamada a la función de residuos por parte de IDAS. No se debe modificar
-- IDAS propone unos valores de las variables y de las derivadas y pide calcular los residuos
     METHOD NO_TYPE funcionResiduos( IN REAL tiempo, IN REAL var[], IN REAL der[], IN INTEGER numPar, IN STRING nomPar[], IN REAL par[], OUT REAL res[] )
			DECLS
			OBJECTS
			BODY
				FOR( i IN 1,numPar )
					setValueReal( nomPar[i], par[i] )
				END FOR
				
				FRES_ARGS( tiempo, var, der, res, FALSE )

		END METHOD

-- llamada a las funciones de cuadraturas por parte de IDAS, para el cálculo de sensibilidades
-- utilizar el vector F_optim que se actualiza al llamar a los residuos
		METHOD NO_TYPE funcionQuadraturas( IN REAL tiempo, IN REAL var[], IN REAL der[], IN INTEGER numPar, IN STRING nomPar[], IN REAL par[], OUT REAL quad[] )
			DECLS
			OBJECTS
			BODY
				FOR( i IN 1,numPar )
					setValueReal( nomPar[i], par[i] )
				END FOR
				
				FRES_ARGS( tiempo, var, der, der, FALSE )
		
		-- funciones de cuadraturas
				FOR( i IN 1,numC +1 )
					quad[i] = F_optim[i]
				END FOR

		END METHOD

--FUNC_PTR<ptrFunRes> ptrRes = optimiz.funcionResiduos
--FUNC_PTR<ptrFunRes> ptrQuad = optimiz.funcionQuadraturas


		METHOD NO_TYPE iniciar()
			DECLS
			OBJECTS
				VECTOR_STRING nombresX
			BODY

				FOR (i IN 1,Nu)
					calculoSens.addU( "uq["+integerToString(i)+"]", uq[i], FALSE)
					calculoSens.addU( "uFr["+integerToString(i)+"]", uFr[i], FALSE)
				END FOR

				SET_INIT_ACTIVE(TRUE)
				EXEC_INIT()
				
				nombresX = getVariablesWithCategory( T_JACOBIAN_VARS )		
				FOR( i IN 1,nombresX.size() )
					calculoSens.addX( nombresX.at(i), getValueReal( nombresX.at(i) ), ( getVarCategory (nombresX.at(i)) == 2 ) )
				END FOR
				numX = nombresX.size()
				
				calculoSens.iniciarRes(TIME, funcionResiduos)
				calculoSens.iniciarQuad(numC + 1, funcionQuadraturas)
				
				calculoSens.setTol( 1e-06 )
				
		END METHOD


     METHOD INTEGER coste_y_restricciones (IN REAL esnopt_x[], IN INTEGER needF, OUT REAL esnopt_F[], IN INTEGER explicit_derivatives, IN INTEGER needG, OUT REAL esnopt_G[])

			DECLS
				REAL arr_quad[ numC + 1 ]						-- array de integrales temporales
				REAL arr_quadsen[ (numC + 1) * numU ]		-- array de sensibilidad de las integrales
				REAL arr_quadsen_p[ (numC + 1) * numU ]	-- array de la derivada temporal de la sensibilidad de las integrales
				
			BODY
				FOR( j IN 1,numC + 1 )
					arr_quad[ j ] = 0
				END FOR
				FOR( j IN 1, (numC + 1) * numU )
					arr_quadsen[ j ] = 0
					arr_quadsen_p[ j ] = 0
				END FOR
				
				FOR (i IN 1,Nu)
					calculoSens.setU( "uq["+integerToString(i)+"]",esnopt_x[MV*(i-1)+1])
					calculoSens.setU( "uFr["+integerToString(i)+"]",esnopt_x[MV*(i-1)+2])
				END FOR

				REL_ERROR = 1e-06 -- transient solver relative tolerance
      		ABS_ERROR = 1e-06 -- transient solver absolute tolerance

				SET_INIT_ACTIVE(TRUE)
					
				TSTOP = Pred_h
				TIME = 0
							
		-- Reiniciar el cálculo de sensibilidades con los nuevos valores de los parámetros				
				calculoSens.reiniciar()
				calculoSens.calculateSens( TSTOP )
				calculoSens.getQuad( arr_quad, arr_quadsen )		-- cálculo de las cuadraturas y sus sensibilidades
				calculoSens.getQuadN( arr_quadsen_p, 1 )			-- cálculo de la derivada de las sensibilidades de las cuadraturas
	
		--	Introduccion de los gradientes de la funcion de costo y restricciones a SNOPT
		--	si la funcion es de camino, hay que introducir la derivada parcial de la cuadratura ( arr_quadsen )
		-- en caso contrario, si la restricción es de punto final, al usar arr_quadsen_p se usa la derivada parcial del valor de la función.

				IF ( explicit_derivatives == 1) THEN
					FOR (i IN 1,numU)
						esnopt_G[i] =arr_quadsen[i]     --arr_quadsen[i]
					END FOR
					FOR (i IN numU+1,(numC+1)*numU)
						esnopt_G[i] = arr_quadsen[i]
					END FOR
					
				END IF

		-- Introducción de los valores de las funcines a SNOPT
		--	si la funcion es de camino, hay que introducir la cuadratura ( arr_quad )
		-- en caso contrario, si la restricción es de punto final, al usar F_optim se usa el valor de la función.
				FOR (i IN 1,numC+1)
					esnopt_F[i] = arr_quad[i]
				END FOR
				
			RETURN 0
			
		END METHOD

		METHOD NO_TYPE control (IN REAL t_sam, IN REAL t_pred, IN REAL state[], IN REAL aux[], IN REAL v_new[], IN REAL error_y[], IN REAL Lu_i[], IN REAL Lu_s[], IN REAL Ly_i[], IN REAL Ly_s[],
	                         IN REAL config[], OUT REAL con[])
			
			  DECLS
	
					REAL dec_var[numU]      	-- valor inicial de las variables de decision, size numU. Necesario para las llamadas a ESnopt
		
	  				REAL xlow[numU]				-- valor inferior de las variables de decision, size numU. Necesario para las llamadas a ESnopt
   				REAL xupp[numU]				-- valor superior de las variables de decision, size numU. Necesario para las llamadas a ESnopt
  
   				REAL Flow[numC + 1]			-- valor inferior de la funcion objetivo y las restricciones, size (numC + función de coste). Necesario para las llamadas a ESnopt
   				REAL Fupp[numC + 1]			-- valor superior de la funcion objetivo y las restricciones, size (numC + función de coste). Necesario para las llamadas a ESnopt
					
					INTEGER calcularSens = 1
					INTEGER infoESnopt = 0
	 	
           OBJECTS
	
		         VECTOR_STRING nombresX		-- variable auxiliar para la inicialización de las sensibilidades
		
				BODY

					t_Sample = t_sam
					Pred_h = t_pred
					
--  Decision variables for this sampling period are initialized to the ones of the 
-- previous period desplaced one sampling time. The index on the 

			      FOR ( i IN 1,numU - MV )  
	    				dec_var[i]= con[MV+i]
					END FOR			
			      FOR ( i IN numU - MV + 1, numU )  
	    				dec_var[i]= con[i]
					END FOR	
					
		-- Inicializaciones de limites de u

				FOR ( j IN 1, MV)
					FOR ( i IN 1,Nu )  
    					xlow[MV*(i-1) + j] = Lu_i[j]    -- limites  de las MV
			      	xupp[MV*(i-1) + j] = Lu_s[j]
					END FOR
				END FOR


				FOR (i IN 1,Nu)
					uq[i] = dec_var[MV*(i-1)+1]
					uFr[i] = dec_var[MV*(i-1)+2]
				END FOR
					
					Ca_inic = state[1]
					Cb_inic = state[2]
					T_inic = state[3]
					Tc_inic = state[4]
					
					uqant = aux[1]
					uFrant = aux[2]
					T0 = aux[3]
					Tc0 = aux[4]	

			      FOR ( i IN 1,Nx )  
	    				v[i]= v_new[i]
					END FOR			
					
					FOR (i IN 1,2)
						error[i] = error_y[i]
					END FOR		
					
				WRITE ("\n\n\t\t------------------- inicializacion en control \n")
				WRITE ("\t\t:\tCa %g \tCb %g \tT %g \tNu %d \tPred_h %g \tt_pred %g\n", Ca_inic, Cb_inic, T_inic, Nu, Pred_h, t_pred)	

				
-- se pasan límites al modelo de las restricciones de camino		

		      LiminfT = Ly_i[1]
  				LimsupT = Ly_s[1]
 				LiminfCb = Ly_i[2]
 				LimsupCb = Ly_s[2]
				
				T_sp = config[1]
				Cb_sp = config[2]
				gamma[1] = config[3]
				gamma[2] = config[4]
				beta[1] = config[5]
				beta[2] = config[6]						

				SET_INIT_ACTIVE(TRUE)
				EXEC_INIT()	
				SET_INIT_ACTIVE(TRUE)
					
				nombresX = getVariablesWithCategory( T_JACOBIAN_VARS )		
				FOR( i IN 1,nombresX.size() )
				calculoSens.setX( nombresX.at(i), getValueReal( nombresX.at(i) ) )
				END FOR
		     		
      	TSTOP = t_pred  --Pred_h
      	CINT = t_sam  --t_Sample
		   TIME = 0        --  Tiempo del modelo
		
		-- Inicializacion del algoritmo de integracion y cálculo de sensibilidades
		-- Inicializacion del algoritmos de optimización
		-- SINTAXIS:
		-- iniciarRes( tiempo inicial desde el que calcular las sensibilidades, puntero a la función de residuos definida previamente )
		-- iniciarQuad( número de cuadraturas a calcular, puntero al cálculo de las funciones de cuadratura )
			calculoSens.iniciarRes(TIME, funcionResiduos)    
			calculoSens.iniciarQuad(numC + 1, funcionQuadraturas)  
		   calculoSens.setTol( 1e-06 )

		-- inicialización de los límites de las restricciones y la función de coste
			Flow[1] = -1.0e10
			Fupp[1] =  1.0e10	
		
		FOR ( i IN 2,numC+1 )  
    			Flow[i] = 0
				Fupp[i] = 0								
		END FOR		

   	--Optimization extern routine call
			setSilentMode(TRUE)
			SET_REPORT_ACTIVE("#MONITOR",FALSE)
		
			esnopt_init (numU, numC)
			esnopt_set_variables_bounds_and_initial_values (xlow, xupp, dec_var)
			esnopt_set_constraints_bounds_and_initial_values (Flow, Fupp, F_optim)
			esnopt_set_cost_function_and_constraints (coste_y_restricciones)
			esnopt_set_explicit_derivatives ( calcularSens )
			esnopt_set_function_precision ( 1.0e-6 )
			esnopt_set_iterations_limit (500) 
			infoESnopt = esnopt ()

		-- Final de la optimización, obtención de los resultados para la simulación.
			setSilentMode(FALSE)
			SET_REPORT_ACTIVE("#MONITOR",TRUE)
		
			esnopt_print_data ()
			esnopt_get_variables_values(dec_var)
			esnopt_free ()	
		
		-- Llamada al integrador
--		RESET()
--			SET_INIT_ACTIVE(TRUE)
		-- Modificación de los parámetros con los nuevos valores obtenidos		
			FOR ( i IN 1,numU )
    			con[i]  = dec_var[i]
			END FOR
			
		END METHOD
		
END CLASS

