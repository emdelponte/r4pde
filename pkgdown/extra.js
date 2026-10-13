// Make dark the default theme. A visitor's own choice in the light switch is
// kept (pkgdown stores it in localStorage under "theme"); only first visits,
// with nothing stored, start in dark mode instead of following the system.
(function () {
  try {
    if (!localStorage.getItem("theme")) {
      localStorage.setItem("theme", "dark");
    }
    if (localStorage.getItem("theme") === "dark") {
      document.documentElement.setAttribute("data-bs-theme", "dark");
    }
  } catch (e) {
    document.documentElement.setAttribute("data-bs-theme", "dark");
  }
})();
