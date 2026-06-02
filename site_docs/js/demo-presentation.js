const initDemoPresentation = () => {
  const button = document.querySelector("[data-demo-presentation-toggle]");

  if (!button) {
    return;
  }

  const activeClass = "demo-presentation-mode";
  const label = button.querySelector("[data-demo-presentation-label]");

  const setActive = (active) => {
    document.documentElement.classList.toggle(activeClass, active);
    document.body.classList.toggle(activeClass, active);
    button.setAttribute("aria-pressed", String(active));

    if (label) {
      label.textContent = active ? "Exit presentation" : "Presentation mode";
    }
  };

  const setPresentationUrl = (active) => {
    const url = new URL(window.location.href);

    if (active) {
      url.searchParams.set("presentation", "1");
    } else {
      url.searchParams.delete("presentation");
    }

    window.history.replaceState({}, "", url);
  };

  const enterPresentation = async () => {
    setActive(true);
    setPresentationUrl(true);

    if (document.fullscreenEnabled && !document.fullscreenElement) {
      try {
        await document.documentElement.requestFullscreen();
      } catch {
        button.dataset.fullscreenBlocked = "true";
      }
    }
  };

  const exitPresentation = async () => {
    if (document.fullscreenElement) {
      await document.exitFullscreen();
    }

    setActive(false);
    setPresentationUrl(false);
  };

  button.addEventListener("click", () => {
    if (document.body.classList.contains(activeClass)) {
      exitPresentation();
    } else {
      enterPresentation();
    }
  });

  document.addEventListener("fullscreenchange", () => {
    if (!document.fullscreenElement && document.body.classList.contains(activeClass)) {
      setActive(false);
      setPresentationUrl(false);
    }
  });

  if (new URLSearchParams(window.location.search).get("presentation") === "1") {
    setActive(true);
  }
};

if (document.readyState === "loading") {
  document.addEventListener("DOMContentLoaded", initDemoPresentation);
} else {
  initDemoPresentation();
}
