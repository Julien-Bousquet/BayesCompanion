#' @import shiny
#' @import graphics
#' @import stats
NULL

#' Compute parameters of beta distribution
#'
#' Compute $a$ and $b$ params of beta distribution
#' from other meaningfull params.
#'
#' @author Julien Bousquet (2022)
#' @param a shape parameter
#' @param b shape parameter
#' @param mean mean of the distribution
#' @param sd standard deviation of the distribution
#' @param k concentration of the distribution
#' @param mode of the distribution
#'
#' @return A list of all possible values of parameters
#' 
#' @examples 
#' betaParams(mean=0.5, k=8)
#' betaParams(mode=0.6, k=8)
#' betaParams(mean=0.5, s=0.2)
#' betaParams(a=12, b=12)
#' 
#' @export
betaParams <- function(a=NA, b=NA, mean=NA, sd=NA, k=NA, mode=NA, verbose=TRUE){
	m <- mean 
	s <- sd
	debg <- FALSE # debugage
	# Everything from a and b
  if(!is.na(a) & !is.na(b)){
		if(debg)cat("a= ",a,'\n')
		if(debg)cat("b= ",b,'\n')
		m <- a/(a+b) # (re-)compute mean from a and b
		if(debg)cat("m= ",m,'\n')
		if(debg)cat("la racine recoit :",m*(1-m)/(a/m+1),'\n')
		res <- list(
				  a=a, b=b, mean=m, mode=(a-1)/(a+b-2), concentration=a+b, 
				  sd=sqrt(m*(1-m)/(a/m+1))
				)
		if(debg)print('List of response is created')
		return(res); #break
  }
	# a and b from mean et sd 
  if(!is.na(m) & !is.na(s)){
	if(verbose==TRUE)print('a and b from mean and sd')
	a <- m*(m*(1-m)/s^2-1)
	b <- (1-m)*(m*(1-m)/s^2-1)
	return(betaParams(a=a, b=b))
  }

	# a and b from mode and concentration
  if(!is.na(mode) & !is.na(k)){
	if(verbose==TRUE)print('a and b from mode and concentration')
	if(k<=2)stop('Concentration of beta distribution must be greater than 2')
	a <- mode*(k-2)+1
	b <- (1-mode)*(k-2)+1
	return(betaParams(a=a, b=b))
  }

	# a and b from mean and concentration
  if(!is.na(m) & !is.na(k)){
	if(verbose==TRUE)print('a and b from mean and concentration')
	a <- m*k
	b <- (1-m)*k
	return(betaParams(a=a,b=b))
  }
}




##########""
# Gamma
#######
#' Compute params of gamma distribution
#'
#' Compute params of gamma distribution
#' from other meaningfull params.
#' Only the couples (shape, rate) or
#' (mean, sd) or (mean, mode) are implemented .
#'
#' @author Julien Bousquet (2022)
#' @param s shape parameter
#' @param r rate parameter
#' @param m mean of the distribution
#' @param s standard deviation of the distribution
#' @param mode of the distribution
#' 
#' @export
#' @examples 
#' gammaParams(m=0.5, mode=0.5) # warning
#' gammaParams(shape=0.6, rate=1) # OK 
#' gammaParams(m=0.5, sd=0.2)     # OK 
#' gammaParams(mode=0.5, sd=5)    # OK
 gammaParams <- function(shape=NA, rate=NA, m=NA, sd=NA, mode=NA, verbose=TRUE){
	debg <- FALSE # debugage
	# Everything from s and r
  if(!is.na(shape) & !is.na(rate)){
		if(shape<=0 | rate<=0)stop("Shape and rate parameter must be positive  for a gamma distribution.")
		if(debg)cat("s= ",shape,'\n')
		if(debg)cat("r= ",rate,'\n')
		res <- list(
			shape=shape, rate=rate, 
			mean=shape/rate, 
			mode=max(c((shape-1)/rate, 0)),
			sd=sqrt(shape)/rate
			)
		if(debg)print('List of response is created')
		return(res); #break
  }
	# shape and rate from mean et sd 
  if(!is.na(m) & !is.na(sd)){
	if(m<=0)stop("The mean parameter m= must be positive for a gamma distribution.")
	if(verbose==TRUE)print('shape and rate from mean and sd')
	shape <- m^2/sd^2
	rate <- m/sd^2
	if(debg)cat(' Debug from mean and sd : shape=', shape, ', rate=', rate,'\n')
	return(gammaParams(shape=shape, rate=rate))
  }

	# shape and rate from mode and sd
  if(!is.na(mode) & !is.na(sd)){
	if(mode <=0)stop("The mode parameter mode= must be positive for a gamma distribution.")
	if(verbose==TRUE)print('shape and rate from mode and sd')
	rate <- (mode+sqrt(mode**2 + 4*sd**2))/2/sd^2
	shape <- 1+mode*rate
	if(debg)cat('rate=', rate,' shape=', shape,'\n')
	return(gammaParams(shape=shape, rate=rate))
  }

  warning("This calculus is not yet implemented.\n")
}




###
# Shiny Applications
#####
# Normale
serverNormale <- function(input, output) {

  mu <- shiny::reactive({return(input$mu)})
  sigma <- shiny::reactive({return(input$sigma)   })


# le graphique de la normale mu, sigma
  output$density <- shiny::renderPlot({
	x <- seq(-10,10, length.out=100)
	plot(x,dnorm(x,mu(),sigma()), type='l', lwd=5, col='blue',
		xlab='x', ylab="p(x)",
		ylim=c(0,1.2),
		main=bquote(atop("Densite de la loi normale", "N(" * mu *'='*.(round(mu(),1)) * ', ' * sigma *'=' *.(round(sigma(),1))  * ')' ) )
	)
	lines(rep(mu()-2*sigma(),2), c(dnorm(mu()-2*sigma(),mu(),sigma())*3,0), 
		col='blue', lwd=3, lty='dotted')
	lines(rep(mu()+2*sigma(),2), c(dnorm(mu()-2*sigma(),mu(),sigma())*3,0), 
		col='blue', lwd=3, lty='dotted')
	text( mu()-sigma(),0.05, "95%", col='blue', cex=1)
	arrows(x0=mu()-2*sigma(),y0=0, x1=mu()+2*sigma(),y1=0, col='blue', code=3, angle=17.5, length=0.15)
	lines(c(mu(),mu()), c(0,1), col='blue', lwd=3)
	text(mu(), 1.05, bquote(mu), col='blue', cex=2)

    })

  output$precision <- shiny::renderPlot({	
	x <- seq(0.1, 5, length.out=100)
	y <- 1/x^2
	preci <- 1/sigma()^2
	plot(x,y, type='l', col='blue', lwd=2,
		xlab=bquote(sigma), ylab='Precision', 
		main=bquote('Precision ='~frac(1,sigma^2 )),
		log='y',
		xlim=c(0,6)
	)
	points(sigma(), preci, col='red', cex=4, pch=16)
	lines(c(0,sigma()), rep(preci,2), lty='dotted',col='blue')
	text(sigma()+0.5, preci, round(preci,3), cex=1.5)
  })
}

uiNormale <- shiny::fluidPage(
# App title ----
  shiny::titlePanel(shiny::HTML("Repr&eacute;sentation de la densit&eacute; de la loi normale")),

  ## Sidebar layout with input and output definitions ----
  shiny::sidebarLayout(
    # Sidebar panel for inputs ----
    shiny::sidebarPanel(
      shiny::sliderInput(inputId = "mu",
                  label = shiny::HTML("&mu;, l'esp&eacute;rance  <br/> param&egrave;tre de position:"),
                  min   = -10,
                  max   = 10,
                  value = 0,
	  step   = 0.1 
	),
      shiny::sliderInput(inputId = "sigma",
                  label =  shiny::HTML("&sigma;, l'&eacute;cart type <br/> param&egrave;tre de dispersion:"),
                  min   = 0.0,
                  max   = 5,
                  value = 1,
	  step   = 0.1 
	),

        shiny::plotOutput(outputId="precision")
   ), 
    shiny::mainPanel(
      # Output: Density ----
      shiny::plotOutput(outputId = "density"),
    )
  )
)

########### server and ui of Beta distribution #########
serverBeta <- function(input, output) {

  mu <- shiny::reactive({return(input$mu)})
  sigma <- shiny::reactive({return(input$sigma)   })
  a <- shiny::reactive({return(ifelse(input$a>0, input$a, 0.001))})
  b <- shiny::reactive({return(ifelse(input$b>0, input$b, 0.001))})
  kappa <- shiny::reactive({return(ifelse(input$kappa>2, input$kappa,2.001))})
  omega <- shiny::reactive({return(input$omega)})


#  graphic of beta mu, sigma
  output$density2 <- shiny::renderPlot({
  	# from mu and sigma
	mu2 <- mu()
	mu2 <- mu()
	sigma2 <- sigma()
	a2 <- mu()*(mu()*(1-mu())/sigma()^2-1)
	b2 <- (1-mu())*(mu()*(1-mu())/sigma()^2-1)
	kappa2 <- a2+b2
	omega2 <- -1
	if( a2 > 1 & b2 > 1) omega2 <- (a2-1)/(a2+b2-2)
	
	x <- seq(0.01,0.99, length.out=100)
	plot(x,dbeta(x,a2,b2), type='l', lwd=5, col='blue',
		xlab='x', ylab="p(x)",
		ylim=c(0,10),
		main=bquote(atop("Densite de la loi beta "*mu==.(mu2)*', '* sigma*'='*.(sigma2), "beta( a ="*.(round(a2,1)) * ', b = '*.(round(b2,1))  * ')' ) )
	)
	aa <- qbeta(0.025,a2,b2) ; bb <- qbeta(0.975,a2,b2)
	lines(rep(aa,2), c(dbeta(aa,a2,b2)+0.5,0), 
		col='blue', lwd=3, lty='dotted')
	lines(rep(bb,2), c(dbeta(bb,a2,b2)+0.5,0), 
		col='blue', lwd=3, lty='dotted')
	text( mu()-sigma(),0.5, "95%", col='blue', cex=1, adj=0)
	arrows(x0=aa,y0=0, x1=bb,y1=0, col='blue', code=3, angle=17.5, length=0.15)
	lines(c(mu(),mu()), c(0,9), col='blue', lwd=3)
	text(mu(), 9.2, bquote(mu), col='blue', cex=2, pos=ifelse(mu2<omega2,2,4))

	# Indication du mode
	if(omega2 <=1 & omega2 >=0){
		lines(c(omega2,omega2), c(0,9), col='red', lwd=3)
		text(omega2, 9.2, bquote(omega), col='red', cex=2, pos=ifelse(mu2<omega2,4,2))
	}
    })

# graphic of beta a ,b
  output$density3 <- shiny::renderPlot({
  	# from a and b
  	a3 <- a()
	b3 <- b()
	kappa3 <- a3+b3
	mu3 <- a3/kappa3
	sigma3 <- sqrt(mu3*(1-mu3)/(a3 + b3 + 1))
	omega3 <- -1
	if( a3>1 & b3>1) omega3 <- (a3-1)/(a3+b3-2)
	
	x <- seq(0.01,0.99, length.out=100)
	plot(x,dbeta(x,a3,b3), type='l', lwd=5, col='blue',
		xlab='x', ylab="p(x)",
		ylim=c(0,10),
		main=bquote(atop("Densite de la loi beta", "beta( a ="*.(round(a3,1)) * ', b = '*.(round(b3,1))  * ')' ) )
	)

	# IDC a 95%
	aa <- qbeta(0.025,a3,b3) ; bb <- qbeta(0.975,a3,b3)
	lines(rep(aa,2), c(dbeta(aa,a3,b3)+0.5,0), 
		col='blue', lwd=3, lty='dotted')
	lines(rep(bb,2), c(dbeta(bb,a3,b3)+0.5,0), 
		col='blue', lwd=3, lty='dotted')
	text( mu3-sigma3,0.5, "95%", col='blue', cex=1, adj=0)
	arrows(x0=aa,y0=0, x1=bb,y1=0, col='blue', code=3, angle=17.5, length=0.15)

	# Indication de l'esperance
	lines(c(mu3,mu3), c(0,9), col='blue', lwd=3)
	text(mu3, 9.2, bquote(mu), col='blue', cex=2, pos=ifelse(mu3<omega3,2,4))
	
	# Indication du mode
	if(omega3 <=1 & omega3 >=0){
		lines(c(omega3,omega3), c(0,9), col='red', lwd=3)
		text(omega3, 9.2, bquote(omega), col='red', cex=2, pos=ifelse(mu3<omega3,4,2))
	}
    })

# graphic of beta kappa, omega
  output$density4 <- shiny::renderPlot({
	kappa4 <- kappa()
	omega4 <- omega()
	a4 <- omega4*(kappa4-2)+1
	b4 <- (1-omega4)*(kappa4-2)+1
	mu4 <- a4/kappa4
	sigma4 <-  sqrt(mu4*(1-mu4) / (a4/mu4 +1)) 

	x <- seq(0.01,0.99, length.out=100)
	plot(x,dbeta(x,a4,b4), type='l', lwd=5, col='blue',
		xlab='x', ylab="p(x)",
		ylim=c(0,10),
		main=bquote(atop("Densite de la loi beta", "beta( a ="*.(round(a4,1)) * ', b = '*.(round(b4,1))  * ')' ) )
	)

	# IDC a 95%
	aa <- qbeta(0.025,a4,b4) ; bb <- qbeta(0.975,a4,b4)
	lines(rep(aa,2), c(dbeta(aa,a4,b4)+0.5,0), 
		col='blue', lwd=3, lty='dotted')
	lines(rep(bb,2), c(dbeta(bb,a4,b4)+0.5,0), 
		col='blue', lwd=3, lty='dotted')
	text( mu4-sigma4,0.5, "95%", col='blue', cex=1, adj=0)
	arrows(x0=aa,y0=0, x1=bb,y1=0, col='blue', code=3, angle=17.5, length=0.15)
	
	
	# Indication de l'esperance
	lines(c(mu4,mu4), c(0,9), col='blue', lwd=3)
	text(mu4, 9.2, bquote(mu), col='blue', cex=2, pos=ifelse(mu4<omega4,2,4))
	
	# Indication du mode
	lines(c(omega4,omega4), c(0,9), col='red', lwd=3)
	text(omega4, 9.2, bquote(omega), col='red', cex=2, pos=ifelse(mu4<omega4,4,2))
	
    })


  output$precision <- shiny::renderPlot({	
	x <- seq(0.1, 5, length.out=100)
	y <- 1/x^2
	preci <- 1/sigma()^2
	plot(x,y, type='l', col='blue', lwd=2,
		xlab=bquote(sigma), ylab='Precision', 
		main=bquote('Precision ='~frac(1,sigma^2 )),
		log='y',
		xlim=c(0,6)
	)
	points(sigma(), preci, col='red', cex=4, pch=16)
	lines(c(0,sigma()), rep(preci,2), lty='dotted',col='blue')
	text(sigma()+0.5, preci, round(preci,3), cex=1.5)
  })
}

uiBeta <- shiny::fluidPage(
# App title ----
  shiny::titlePanel(shiny::HTML("Repr&eacute;sentation de la densit&eacute; de la loi b&ecirc;ta :")),
  #plotOutput(outputId = "density2"),
	  
  shiny::fluidRow(
  ############################ col gauche mu et sigma
    shiny::column(4,
      shiny::sliderInput(inputId = "mu",
                  label = shiny::HTML("&mu;, l'esp&eacute;rance, param&egrave;tre de position:"),
                  min   = 0,
                  max   = 1,
                  value = 0.5,
	  step   = 0.05 
	),
	shiny::br(),	

      shiny::sliderInput(inputId = "sigma",
                  label =  shiny::HTML("&sigma;, l'&eacute;cart type, param&egrave;tre de dispersion:"),
                  min   = 0.0,
                  max   = 0.5,
                  value = 0.25,
	  step   = 0.05 
	),
	shiny::br(),
	
      # Output: Density ----
      shiny::plotOutput(outputId ="density2", width = "60%")
      ), # premiere colonne
 
 
 ############################# colonne centrale a et b
    shiny::column(4, 
	shiny::sliderInput(inputId = "a",
                  label = shiny::HTML("a, param&egrave;tre a :"),
                  min   = 0.0,
                  max   = 50,
                  value = 1,
	  step   = 0.1 
	),
	shiny::br(),
	
      shiny::sliderInput(inputId = "b",
                  label =  shiny::HTML("b, param&egrave;tre b :"),
                  min   = 0.0,
                  max   = 50,
                  value = 1,
	  step   = 0.1 
	),
	shiny::br(),
	
      shiny::plotOutput(outputId = "density3", width = "60%")
      ), # colonne centrale a et b


############################# colonne droite kappa, omega
	shiny::column(4,
	shiny::sliderInput(inputId = "omega",
                  label =  shiny::HTML("&omega;, le mode, param&egrave;tre de position:"),
                  min   = 0,
                  max   = 1,
                  value = 0.5,
	  step   = 0.1 
	),
	shiny::br(),

	shiny::sliderInput(inputId = "kappa",
                  label =  shiny::HTML("&kappa;, la concentration , param&egrave;tre de dispersion :"),
                  min   = 2,
                  max   = 100,
                  value = 2,
	  step   = 1 
	),
	shiny::br(),
	
      shiny::plotOutput(outputId = "density4", width = "60%")
   ) 

    
  )
)


##### server and ui of student's distribution #########
serverStudent <- function(input, output) {
  nu <- shiny::reactive({return(input$nu)})  
  mu <- shiny::reactive({return(input$mu)})
  sigma <- shiny::reactive({return(input$sigma)   })


# le graphique de la normale mu, sigma
  output$density <- shiny::renderPlot({
	x <- sort(c(seq(-10,10, length.out=100), 
	            seq(mu()-sigma(), mu()+sigma(), length.out=100)))
	y <- ggdist::dstudent_t(x,df=nu(), mu=mu(), sigma=sigma())
	plot(x,dnorm(x, mean=mu(), sd=sigma()), type='l', lwd=2, col='grey',
		xlab='x', ylab="p(x)",
		ylim=c(0,1.2),
		main=bquote(atop("Densite de la loi de Student", 
		  "t("*nu*'='*.(round(nu(),1))*', '*mu *'='*.(round(mu(),1)) * ', ' * sigma *'=' *.(round(sigma(),1))  * ')' ) )
	)
	
	lines(x,y, type='l', lwd=5, col='blue')
	a <- ggdist::qstudent_t(0.025,df=nu(), mu=mu(), sigma=sigma())
	b <- ggdist::qstudent_t(0.975,df=nu(), mu=mu(), sigma=sigma())
	# IDC de t
	lines(rep(a,2), c(ggdist::dstudent_t(a,df=nu(), mu=mu(), sigma=sigma())+0.15,0), 
		col='blue', lwd=3, lty='dotted')
	lines(rep(b,2), c(ggdist::dstudent_t(b,df=nu(), mu=mu(), sigma=sigma())+0.1,0), 
		col='blue', lwd=3, lty='dotted')
	# IDC de la normale
	lines(rep(mu()-2*sigma(),2), dnorm(mu()-2*sigma(), mu(), sigma())+c(-0.05, +0.05), 
		col='grey', lwd=3, lty='dotted')
	lines(rep(mu()+2*sigma(),2), dnorm(mu()+2*sigma(), mu(), sigma())+c(-0.05, +0.05), 
		col='grey', lwd=3, lty='dotted')

	text( mu()+sigma(),0.05, "95%", col='blue', cex=1, adj=0)
	arrows(x0=a,y0=0, x1=b,y1=0, col='blue', code=3, angle=17.5, length=0.15)
	lines(c(mu(),mu()), c(0,1), col='blue', lwd=3)
	text(mu(), 1.05, bquote(mu), col='blue', cex=2)
	text(mu()+2*sigma(), max(y)*0.86, adj=0, cex=1.5, col='blue',
	    bquote("t("*nu*'='*.(round(nu(),1))*', '*mu *'='*.(round(mu(),1)) * ', ' * sigma *'=' *.(round(sigma(),1))  * ')'  ) )
	text(mu()-2*sigma(), max(y)*0.86, adj=1, cex=1.5, col='grey',
	    bquote("N("*mu *'='*.(round(mu(),1)) * ', ' * sigma *'=' *.(round(sigma(),1))  * ')'  ) )
	    arrows(mu()+2*sigma(), max(y)*0.86, mu()+sigma(), 
	        ggdist::dstudent_t(mu()+sigma(),df=nu(), mu=mu(), sigma=sigma()), 
	        col='blue', angle=17.5, length=0.15)
	    arrows(mu()-2*sigma(), max(y)*0.86, mu()-sigma(), 
	        dnorm(mu()-sigma(), mean=mu(), sd=sigma()), 
	        col='grey', angle=17.5, length=0.15)
	 
    })

  output$precision <- shiny::renderPlot({	
	x <- seq(0.1, 5, length.out=100)
	y <- 1/x^2
	preci <- 1/sigma()^2
	plot(x,y, type='l', col='blue', lwd=2,
		xlab=bquote(sigma), ylab='Precision', 
		main=bquote('Precision ='~frac(1,sigma^2 )),
		log='y',
		xlim=c(0,6)
	)
	points(sigma(), preci, col='red', cex=4, pch=16)
	lines(c(0,sigma()), rep(preci,2), lty='dotted',col='blue')
	text(sigma()+0.5, preci, round(preci,3), cex=1.5)
  })
}

uiStudent <- shiny::fluidPage(
# App title ----
  shiny::titlePanel(shiny::HTML("Repr&eacute;sentation de la densit&eacute; de la loi de Student")),

  # Sidebar layout with input and output definitions ----
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::sliderInput(inputId = "nu",
                  label = shiny::HTML("&nu;, la normalit&eacute;<br/> degr&eacute; de libert&eacute; :"),
                  min   = 1,
                  max   = 30,
                  value = 2,
	  step   = 0.1 
	),

      shiny::sliderInput(inputId = "mu",
                  label =shiny::HTML("&mu;, l'esp&eacute;rance  <br/> param&egrave;tre de position:"),
                  min   = -10,
                  max   = 10,
                  value = 0,
	  step   = 0.1 
	),

      shiny::sliderInput(inputId = "sigma",
                  label =  shiny::HTML("&sigma;, param&egrave;tre de dispersion :<br/>(n'est pas &eacute;gal &agrave; l'&eacute;cart type)"),
                  min   = 0.0,
                  max   = 5,
                  value = 1,
	  step   = 0.1 
	),

        shiny::plotOutput(outputId="precision")
   ), 

    shiny::mainPanel(
      # Output: Density ----
      shiny::plotOutput(outputId = "density"),
    )
  )
)


################################
# Call to shiny functions  ################
###############################
#' Shiny application for normal distribution
#'
#' Demonstrate behaviour of the parameters of 
#' normal distribution, in a shiny application.
#' @return None
#' @author Julien Bousquet (2022)
#' @examples 
#' normaleShiny()
#'
#' @export
normaleShiny <- function(){
	shiny::shinyApp(ui = uiNormale, server = serverNormale)
}
#' Demonstrate behaviour of the parameters of 
#' beta distribution, in a shiny application.
#'
#' @author Julien Bousquet (2022)
#' @return None
#' @examples 
#' betaShiny()
#'
#' @export
betaShiny <-  function(){
	shiny::shinyApp(ui = uiBeta, server = serverBeta)
}
#' Demonstrate behaviour of the parameters of 
#' Student distribution, in a shiny application.
#'
#' @author Julien Bousquet (2022)
#' @return None
#' @examples 
#' StudentShiny()
#'
#' @export

StudentShiny <-  function(){
	shiny::shinyApp(ui = uiStudent, server = serverStudent)
}