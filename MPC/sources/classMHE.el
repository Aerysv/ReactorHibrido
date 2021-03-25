/*-----------------------------------------------------------------------------------------
 LIBRARY: PLANTAHIBRIDA
 FILE: classMHE
 CREATION DATE: 17/07/2020
-----------------------------------------------------------------------------------------*/
USE MATH
USE OPTIM_METHODS
		
CLASS estimMHE IS_A CostoMHE_par1

	DECLS

		CONST INTEGER numC = 4			-- número de restricciones del problema de optimización
		CONST INTEGER numU = 8			-- número de variables de decisión  MV*Nu + slacks
		INTEGER numX = 4					-- número de variables de estado del modelo	
		REAL t_mhe = 2
		INTEGER MV = 2 					-- numero de variables manipuladas
		
		
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


				SET_INIT_ACTIVE(TRUE)
				EXEC_INIT()	
				
				nombresX = getVariablesWithCategory( T_JACOBIAN_VARS )		
				FOR( i IN 1,nombresX.size() )
					calculoSens.addX( nombresX.at(i), getValueReal( nombresX.at(i) ), ( getVarCategory (nombresX.at(i)) == 2 ) )
				END FOR
				numX = nombresX.size()
				
		-- addU( nombre de la variable, valor inicial de la misma, booleano que indica si el parámetro es un valor inicial )				
				FOR( i IN 1,numX )
					calculoSens.addU("v["+integerToString(i)+"]", v[i], FALSE )
				END FOR

				calculoSens.addU( "Ca_N", Ca_N, TRUE)
				calculoSens.addU( "Cb_N", Cb_N, TRUE)
				calculoSens.addU( "T_N",  T_N, TRUE)				
				calculoSens.addU( "Tc_N", Tc_N, TRUE)		
				
				calculoSens.iniciarRes(TIME, funcionResiduos)
				calculoSens.iniciarQuad(numC + 1, funcionQuadraturas)
				
				calculoSens.setTol( 1e-06 )
				
		END METHOD


     METHOD INTEGER coste_y_restricciones (IN REAL esnopt_x[], IN INTEGER needF, OUT REAL esnopt_F[], IN INTEGER explicit_derivatives, IN INTEGER needG, OUT REAL esnopt_G[])

			DECLS
				REAL arr_quad[ numC + 1 ]						-- array de integrales temporales
				REAL arr_quadsen[ (numC + 1) * numU ]		-- array de sensibilidad de las integrales
				REAL arr_quadsen_p[ (numC + 1) * numU ]	-- array de la derivada temporal de la sensibilidad de las integrales
			OBJECTS
				-- VECTOR_STRING nombresX
				--EIDAS calculoSens
				
			BODY
				FOR( j IN 1,numC + 1 )
					arr_quad[ j ] = 0
				END FOR
				FOR( j IN 1, (numC + 1) * numU )
					arr_quadsen[ j ] = 0
					arr_quadsen_p[ j ] = 0
				END FOR
				

				FOR( i IN 1,numX )
					calculoSens.setU( "v["+integerToString(i)+"]",esnopt_x[i])
				END FOR
				  calculoSens.setU( "Ca_N", esnopt_x[5] )
				  calculoSens.setU( "Cb_N", esnopt_x[6] )
				  calculoSens.setU( "T_N", esnopt_x[7] )
				  calculoSens.setU( "Tc_N", esnopt_x[8] )				  

				REL_ERROR = 1e-06 -- transient solver relative tolerance
      		ABS_ERROR = 1e-06 -- transient solver absolute tolerance
      	--	TOLERANCE = 1e-07 -- steady solver relative tolerance
--				SET_INIT_ACTIVE(TRUE)
--				EXEC_INIT()	
				SET_INIT_ACTIVE(TRUE)
					
--				nombresX = getVariablesWithCategory( T_JACOBIAN_VARS )		
--				FOR( i IN 1,nombresX.size() )
--					calculoSens.setX( nombresX.at(i), getValueReal( nombresX.at(i) ) )
--				END FOR

				TSTOP = t_mhe*Ne
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
						esnopt_G[i] =arr_quadsen[i]     --arr_quadsen_p[i]
					END FOR
					FOR (i IN numU+1,(numC+1)*numU)
						esnopt_G[i] = arr_quadsen[i]
					END FOR
					
				END IF

		-- Introducción de los valores de las funcines a SNOPT
		--	si la funcion es de camino, hay que introducir la cuadratura ( arr_quad )
		-- en caso contrario, si la restricción es de punto final, al usar F_optim se usa el valor de la función.
				esnopt_F[1] = arr_quad[1]		--F_optim[1]
				esnopt_F[2] = arr_quad[2]
				esnopt_F[3] = arr_quad[3]
				esnopt_F[4] = arr_quad[4]
				esnopt_F[5] = arr_quad[5]
				
			RETURN 0
			
		END METHOD

		METHOD NO_TYPE mhe (IN REAL t_sam, IN REAL x_ant[], IN REAL accion[], IN REAL perturbacion[], IN REAL medida[], IN REAL u_ant[], IN REAL margen_v, IN REAL margen_x, IN REAL beta_v, IN REAL beta_N, OUT REAL state[], OUT REAL pert[], OUT REAL error_y[] )
		
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
--		         VECTOR_STRING nombresU		-- variable auxiliar para la inicialización de las sensibilidades
		
				BODY

					t_mhe = t_sam
					beta_xv = beta_v 
					beta_xN = beta_N					
					
--  Decision variables for this sampling period are initialized 

					FOR(i IN 1, numX)
				    	dec_var[i]= pert[i]
						--dec_var[i] = 0.01
						dec_var[numX + i]= x_ant[i]
					END FOR
				
					FOR(i IN 1, numX)
				    	v_ant[i]= pert[i]
					END FOR	
					
		-- Inicializaciones de limites de v y x_N

					FOR (i IN 1, numX)
  						xlow[i] = - margen_v   -- limites  de las v
			      	xupp[i] =   margen_v
						
  						xlow[numX + i] =  x_ant[i]*(1 - margen_x/100)    -- limites  de las x_N
			      	xupp[numX + i] =  x_ant[i]*(1 + margen_x/100)						
					END FOR
					
					Ca_N = x_ant[1]		
					Cb_N = x_ant[2]		
					T_N  = x_ant[3]	
					Tc_N = x_ant[4]		

					Ca_Ne = x_ant[1]		
					Cb_Ne = x_ant[2]		
					T_Ne  = x_ant[3]		
					Tc_Ne = x_ant[4]		

				FOR (j IN 1,Ne+1)
					-- Medidas de las acciones de control
					setValueReal( "tq["+integerToString(j)+"]",accion[MV*(j-1)+1])
					setValueReal( "tFr["+integerToString(j)+"]",accion[MV*(j-1)+2])
					-- Perturbaciones medidas
					setValueReal( "tT0["+integerToString(j)+"]",perturbacion[Nd*(j-1)+1])
					setValueReal( "tTc0["+integerToString(j)+"]",perturbacion[Nd*(j-1)+2])
					-- Estados medidos
					setValueReal( "tCa["+integerToString(j)+"]",medida[Nm*(j-1)+1])
					setValueReal( "tCb["+integerToString(j)+"]",medida[Nm*(j-1)+2])
					setValueReal( "tT["+integerToString(j)+"]",medida[Nm*(j-1)+3])
					setValueReal( "tTc["+integerToString(j)+"]",medida[Nm*(j-1)+4])
				END FOR					
		
				WRITE ("\n\n\t\t------------------- inicializacion en mhe \n")
				WRITE ("\t\t:\tCa %g \tCb %g \tT %g \tTc %g \n", Ca_N, Cb_N, T_N, Tc_N)	

				
-- se pasan límites al modelo de las path constraints		

		-- inicialización de los límites de las restricciones y la función de coste
				Flow[1] = -1.0e6
				Fupp[1] =  1.0e6			
		
				FOR ( i IN 2,numC+1 )  
    				Flow[i] = 0
					Fupp[i] = 0								
				END FOR	
		
				
				SET_INIT_ACTIVE(TRUE)
				EXEC_INIT()	
				SET_INIT_ACTIVE(TRUE)
					
				nombresX = getVariablesWithCategory( T_JACOBIAN_VARS )		
				FOR( i IN 1,nombresX.size() )
				calculoSens.setX( nombresX.at(i), getValueReal( nombresX.at(i) ) )
				END FOR
		
      -- simulate a transient in range[TIME,TSTOP] reporting every CINT
     		
      		TSTOP = t_mhe*Ne   --Pred_h
      		CINT = t_mhe       --t_Sample
		   	TIME = 0           --  Tiempo del modelo
		
		-- Inicializacion del algoritmo de integracion y cálculo de sensibilidades
			calculoSens.iniciarRes(TIME, funcionResiduos)    
			calculoSens.iniciarQuad(numC + 1, funcionQuadraturas) 
		   calculoSens.setTol( 1e-06 )

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
			FOR ( i IN 1,numX )
    			pert[i] = dec_var[i]
			END FOR
					
			FOR ( i IN 1,numX )
					v[i] = pert[i]
			END FOR		
			
			Ca_N = dec_var[5]
			Cb_N = dec_var[6]
			T_N  = dec_var[7]
			Tc_N = dec_var[8]
								
      	TSTOP = t_mhe*Ne  -- hasta t
      	CINT = t_mhe  --t_Sample
		  	TIME = 0        --  Tiempo del modelo
			SET_INIT_ACTIVE(TRUE)
			
			INTEG()
			
			state[1] = Ca
			state[2] = Cb									
			state[3] = T
			state[4] = Tc

			FOR( i IN 1,2 )
			  	error_y[i]= error[i]
			END FOR				
			
		
			WRITE ("\n\n\t\t------------------- salida de mhe \n")
			WRITE ("\t\t:\tCa %g \tCa_m %g \tCa_N %g \tpert[1] %g \tstate[1] %g \n",Ca, Ca_m, Ca_N, pert[1],state[1])
			WRITE ("\t\t:\tCb %g \tCb_m %g \tCb_N %g \tpert[2] %g \tstate[2] %g \n",Cb, Cb_m, Cb_N, pert[2],state[2])	
			WRITE ("\t\t:\tT %g \tT_m %g \tT_N %g \tpert[3] %g \tstate[3] %g \n",T, T_m, T_N, pert[3],state[3])	
			WRITE ("\t\t:\tTc %g \tTc_m %g \tTc_N %g \tpert[4] %g \tstate[4] %g \n",Tc, Tc_m, Tc_N, pert[4],state[4])	
			WRITE ("\t\t:\tJ_costo %g \tJ_costo_m %g \tJ_costo_N %g \tJ_costo_v %g \tbeta_xv %g \n",J_costo, J_costo_m, J_costo_N, J_costo_v, beta_xv)

		END METHOD
		
END CLASS

