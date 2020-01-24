#### compare deltaG method results

#method 1
method1 = fread("dataset/dG/method1_bestmodel.txt")
method1[,id := paste0(Pos,Mut,collapse=""),.(Pos,Mut)]
method1_pars = fread("dataset/dG/method1_globalpar.txt")
#method 2
method2 = fread("dataset/dG/method2_top100vars_bestmodel.txt")
method2_pars = fread("dataset/dG/method2_top100vars_bestpar.txt")


X = merge(method1[,.(id,f_ddG_method1 = f_ddG_local,b_ddG_method1 = b_ddG_local,type)],
      method2[Nmut==1,.(id=id1,f_ddG_method2 = f_ddG_local1,b_ddG_method2 = b_ddG_local1,f_fitness,b_fitness)])

ggpairs(X,columns = grep("[fb]",names(X),value=T),aes(color=type,alpha=0.5))

f_dg=seq(-3,3,0.1)
ggplot(X,aes(f_ddG_method2,f_fitness)) +
  geom_point(aes(color=type,alpha=0.5)) + 
  geom_line(inherit.aes = F,data=data.table(x=seq(-2,3,0.1),
                                            y=function_folding_dG2F(seq(-2,3,0.1),method2_pars$f_dGwt,method2_pars$f_bgr,method2_pars$f_scale)),
            aes(x,y))
  
ggplot(X,aes(f_ddG_method2,b_fitness)) +
  geom_point(aes(color=type,alpha=0.5)) + 
  geom_line(inherit.aes = F,data=data.table(x=seq(-2,3,0.1),
                                            y=function_binding_dG2F(b_dG = 0,f_dG = seq(-2,3,0.1),
                                                                    b_dGwt = method2_pars$b_dGwt,
                                                                    f_dGwt = method2_pars$f_dGwt,
                                                                    b_bgr = method2_pars$b_bgr,
                                                                    b_scale = method2_pars$b_scale)),
            aes(x,y))
