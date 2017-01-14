setGeneric("get_curve",
	def = function( object, par, times, options=list() ){
		standardGeneric("get_curve");
	})

setGeneric("get_param_info",
	def = function(object, times, options=list()){
		standardGeneric("get_param_info");
	})


setGeneric("check_param",
	def = function(object, par, times, options=list()){
		standardGeneric("check_param");
	})

setGeneric("get_simu_param",
	def = function(object, times, options=list()){
		standardGeneric("get_simu_param");
	})

setGeneric("est_init_param",
	def = function( object, pheY, pheX, pheT, options=list()){
		standardGeneric("est_init_param");
	})


setGeneric("get_matrix",
	def = function( object, par, times, options=list() ){
		standardGeneric("get_matrix");
	})
