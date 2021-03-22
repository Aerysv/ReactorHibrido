/*-----------------------------------------------------------------------------------------
 LIBRARY: MATH_ISA
 FILE: SIGMOIDE
 AUTHOR: DPTO. INGENIERÍA DE SISTEMAS Y AUTOMÁTICA
 COMPANY: UNIVERSIDAD DE VALLADOLID
 CREATION DATE: 09/03/2017
-----------------------------------------------------------------------------------------*/

FUNCTION REAL sigmoide( REAL EQUIS, REAL CENTRO, REAL SIGMA)
//	DECLS
//		REAL Z
	BODY
//		Z = (1 / (1 + exp(-(EQUIS- CENTRO) * SIGMA)))
//		RETURN Z
		RETURN (1 / (1 + exp(-(EQUIS- CENTRO)/60 * SIGMA)))
END FUNCTION
