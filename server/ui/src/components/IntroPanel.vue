<script setup>
import ExampleImage from '../assets/image.png'
</script>

<template>
  <div class="h-full grid">
    <div class="col-5 flex-center">
      <div>
        <h1 class="mt-0 mb-2 text-6xl">Data Table Plugin</h1>
        <div class="mb-8 text-500 text-xl">
          View your multi-frame molecule metadata<br />
          in an interactive table.
        </div>

        <div
          class="surface-card px-7 py-5 text-center"
          style="border-radius: 8px"
        >
          <slot />
        </div>
      </div>
    </div>

    <div class="col-7 px-6 flex-center surface-50 text-center">
      <div>
        <h2 class="text-4xl">Example Output</h2>
        <img :src="ExampleImage" alt="example" class="w-full" />
        <div class="mt-4 text-500 text-xl">
          Best used for analyzing metadata for a multi-frame SDF small molecule
          entry
        </div>
      </div>
    </div>
  </div>
</template>
