{
  gROOT->Reset();
  TChain chain("h1");
  chain.Add("../Ds2317.root");

  Int_t
    BS_0      = chain.GetEntries("cn_dsii!=2 && cn_dsi!=2"),
    S_0       = chain.GetEntries("cn_dsii!=2 && cn_dsi!=2 && abs(mc_ds17)==10431"),
    BS_pst17  = chain.GetEntries("cn_dsii!=2 && cn_dsi!=2 && pst_ds17>2.79"),
    S_pst17   = chain.GetEntries("cn_dsii!=2 && cn_dsi!=2 && abs(mc_ds17)==10431 && pst_ds17>2.79"),
    BS_r2     = chain.GetEntries("cn_dsii!=2 && cn_dsi!=2 && r2>0.205"),
    S_r2      = chain.GetEntries("cn_dsii!=2 && cn_dsi!=2 && abs(mc_ds17)==10431 && r2>0.205"),
    BS_xp     = chain.GetEntries("cn_dsii!=2 && cn_dsi!=2 && x_p>0.21"),
    S_xp      = chain.GetEntries("cn_dsii!=2 && cn_dsi!=2 && abs(mc_ds17)==10431 && x_p>0.21"),
    BS_pstdsp = chain.GetEntries("cn_dsii!=2 && cn_dsi!=2 && pst_dsp>1.8"),
    S_pstdsp  = chain.GetEntries("cn_dsii!=2 && cn_dsi!=2 && abs(mc_ds17)==10431 && pst_dsp>1.8"),
    BS_hel    = chain.GetEntries("cn_dsii!=2 && cn_dsi!=2 && hel_phi>0.44"),
    S_hel     = chain.GetEntries("cn_dsii!=2 && cn_dsi!=2 && abs(mc_ds17)==10431 && hel_phi>0.44");

  cout << "              N/A:    S+B=" << BS_0      << '\t' << "S=" << S_0      << '\t' << "a=S/sqrt(S+B)=" << S_0/sqrt(BS_0)           << '\n';
  cout << "  p*(Ds2317)>2.79:    S+B=" << BS_pst17  << '\t' << "S=" << S_pst17  << '\t' << "a=S/sqrt(S+B)=" << S_pst17/sqrt(BS_pst17)   << '\n';
  cout << "         r2>0.205:    S+B=" << BS_r2     << '\t' << "S=" << S_r2     << '\t' << "a=S/sqrt(S+B)=" << S_r2/sqrt(BS_r2)         << '\n';
  cout << "          xP>0.21:    S+B=" << BS_xp     << '\t' << "S=" << S_xp     << '\t' << "a=S/sqrt(S+B)=" << S_xp/sqrt(BS_xp)         << '\n';
  cout << "      p*(Dsp)>1.8:    S+B=" << BS_pstdsp << '\t' << "S=" << S_pstdsp << '\t' << "a=S/sqrt(S+B)=" << S_pstdsp/sqrt(BS_pstdsp) << '\n';
  cout << "cos|Theta_H|>0.44:    S+B=" << BS_hel    << '\t' << "S=" << S_hel    << '\t' << "a=S/sqrt(S+B)=" << S_hel/sqrt(BS_hel)       << '\n';
}
