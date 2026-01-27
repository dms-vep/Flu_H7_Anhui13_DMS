import { defineConfig } from "vitepress";

// https://vitepress.dev/reference/site-config
export default defineConfig({
  lang: "en-US",
  title: "Pseudovirus deep mutational scanning of H7 hemagglutinin (A/Anhui/1/2013)",
  description:
    "Data and interactive figures for pseudovirus deep mutational scanning of influenza H7 HA",
  base: "/Flu_H7_Anhui13_DMS/",
  appearance: false,
  themeConfig: {
    // https://vitepress.dev/reference/default-theme-config
    nav: [
      { text: "Home", link: "/" },
      { text: "Appendix", link: "/appendix", target: "_self" },
    ],
    socialLinks: [{ icon: "github", link: "https://github.com/dms-vep/Flu_H7_Anhui13_DMS" }],
    footer: {
      message: '<a href="https://doi.org/10.64898/2026.01.05.697808">Ahn, Yu, et al (2026)</a>',
    },
  },
});
