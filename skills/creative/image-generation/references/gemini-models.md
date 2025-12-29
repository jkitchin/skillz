# Google Gemini Image Generation Models

> **Important (December 2025)**: This documentation uses the `google-genai` SDK. The `google-generativeai`
> package has been deprecated. See [migration guide](https://ai.google.dev/gemini-api/docs/migrate).

## Overview

Google's Gemini models include advanced image generation capabilities through the `generateContent` endpoint with image response modalities. These models can create images from text, edit existing images, and maintain context across multi-turn conversations.

## Available Models

### gemini-2.5-flash-image (Nano Banana)

**Specifications:**
- Maximum Resolution: 1024px
- Speed: Fast (optimized for rapid generation)
- Cost: Lower cost per image
- Best Use Cases: Prototyping, high-volume operations, rapid iteration

**Strengths:**
- Quick generation times
- Good for exploring multiple concepts
- Cost-effective for large batches
- Suitable for web/social media resolutions

**Limitations:**
- Limited to 1024px maximum
- Less detail than Pro model
- May struggle with complex compositions

**Ideal For:**
- Draft concepts and iterations
- Social media graphics
- Web banners and assets
- Quick mockups
- Volume generation

### gemini-3-pro-image-preview (Nano Banana Pro)

**Specifications:**
- Maximum Resolution: 4K (4096px)
- Speed: Slower but higher quality
- Cost: Higher cost per image
- Best Use Cases: Professional assets, final deliverables, complex compositions

**Strengths:**
- Up to 4K resolution for professional use
- Better at following complex instructions
- Superior text rendering in images
- Google Search grounding (access real-time data)
- Multiple reference images (up to 14)
- Enhanced detail and quality

**Limitations:**
- Slower generation time
- Higher API costs
- May be overkill for simple requests

**Ideal For:**
- Final deliverables and publications
- Print materials
- Marketing assets
- Designs with text elements
- Complex multi-element compositions
- Current events visualization (with Search grounding)

## API Usage

### Basic Setup

```python
from google import genai
from google.genai import types

# Initialize client (uses GEMINI_API_KEY or GOOGLE_API_KEY env var automatically)
client = genai.Client()
```

### Text-to-Image Generation

```python
# Simple generation
response = client.models.generate_content(
    model="gemini-2.5-flash-image",
    contents="A serene mountain lake at sunset, photorealistic",
    config=types.GenerateContentConfig(
        response_modalities=["TEXT", "IMAGE"]
    )
)

# Save image
for part in response.parts:
    if part.text is not None:
        print(part.text)
    elif part.inline_data is not None:
        image = part.as_image()
        image.save("output.png")
```

### Configuration Options

```python
response = client.models.generate_content(
    model="gemini-2.5-flash-image",
    contents=prompt,
    config=types.GenerateContentConfig(
        response_modalities=["TEXT", "IMAGE"],
        image_config=types.ImageConfig(
            aspect_ratio="16:9",  # Options: 1:1, 3:4, 4:3, 9:16, 16:9, 21:9
            image_size="2K",      # Options: 1K, 2K, 4K (Pro only)
        )
    )
)
```

**Aspect Ratios:**
- `1:1` - Square (1024x1024) - Social media, logos, icons
- `3:4` - Portrait (768x1024) - Mobile screens, portraits
- `4:3` - Landscape (1024x768) - Presentations, standard screens
- `9:16` - Tall portrait (576x1024) - Stories, vertical video
- `16:9` - Widescreen (1024x576) - Video, cinematic
- `21:9` - Ultra-wide (1024x439) - Panoramic, cinematic

**Image Sizes (Pro only):**
- `1K` - ~1024px (fastest)
- `2K` - ~2048px (balanced)
- `4K` - ~4096px (highest quality)

### Image Editing

```python
from PIL import Image

# Load existing image
input_image = Image.open("photo.jpg")

# Edit with text instructions
response = client.models.generate_content(
    model="gemini-2.5-flash-image",
    contents=[input_image, "Change the background to a starry night sky"],
    config=types.GenerateContentConfig(
        response_modalities=["TEXT", "IMAGE"]
    )
)

# Save edited image
for part in response.parts:
    if part.inline_data is not None:
        image = part.as_image()
        image.save("edited.png")
```

### Multi-turn Refinement

```python
# Start chat for iterative refinement
chat = client.chats.create(
    model="gemini-2.5-flash-image",
    config=types.GenerateContentConfig(
        response_modalities=["TEXT", "IMAGE"]
    )
)

# First generation
response1 = chat.send_message("A futuristic city skyline at night")

# Refinement 1
response2 = chat.send_message("Add more neon lights and flying cars")

# Refinement 2
response3 = chat.send_message("Make the scene more cyberpunk with rain reflections")

# Save each iteration
for i, resp in enumerate([response1, response2, response3], 1):
    for part in resp.parts:
        if part.inline_data is not None:
            image = part.as_image()
            image.save(f"iteration_{i}.png")
```

### Multiple Reference Images (Pro Only)

```python
from PIL import Image

# Load reference images
ref1 = Image.open("style_reference.jpg")
ref2 = Image.open("composition_reference.jpg")
ref3 = Image.open("color_reference.jpg")

# Generate with multiple references
response = client.models.generate_content(
    model="gemini-3-pro-image-preview",
    contents=[
        ref1,
        ref2,
        ref3,
        "Create a portrait combining the artistic style of the first image, "
        "the composition of the second, and the color palette of the third"
    ],
    config=types.GenerateContentConfig(
        response_modalities=["TEXT", "IMAGE"],
        image_config=types.ImageConfig(image_size="4K")
    )
)
```

### Google Search Grounding (Pro Only)

```python
# Use Search grounding for current events/real locations
response = client.models.generate_content(
    model="gemini-3-pro-image-preview",
    contents="Create an aerial view visualization of Tokyo's Shibuya crossing "
             "showing current building layouts and architecture",
    config=types.GenerateContentConfig(
        response_modalities=["TEXT", "IMAGE"],
        tools=[{"google_search": {}}]  # Enable Search grounding
    )
)
```

## Prompt Engineering Best Practices

### Structure for Success

**Formula:** `[Subject] + [Style] + [Details] + [Lighting] + [Mood] + [Technical]`

**Examples:**

**Photorealistic:**
```
Portrait of an elderly craftsman, photorealistic style, weathered hands holding
antique tools, warm natural lighting from window, contemplative mood, shot on
85mm lens, shallow depth of field, high detail
```

**Illustration:**
```
Kawaii robot character, cute illustration style, pastel colors with mint and
peach, rounded friendly shapes, bold black outlines, cel-shaded, cheerful
expression, sticker format with white border
```

**Product Photography:**
```
Premium coffee maker, professional product photography, brushed stainless steel
finish, white seamless background, studio lighting with soft shadows, 3/4 angle
view, clean minimal composition, commercial quality
```

**Logo Design:**
```
Tech company logo, minimalist geometric design, interconnected nodes forming
abstract network, blue gradient with silver accents, professional and modern,
vector style, clean lines, suitable for digital and print
```

### Effective Techniques

1. **Be Specific About Style:**
   - Good: "photorealistic portrait, 85mm lens, golden hour lighting"
   - Bad: "nice portrait"

2. **Include Technical Details:**
   - Camera: "shot on 24mm wide angle", "macro photography", "drone aerial view"
   - Lighting: "rim light", "soft diffused", "dramatic side lighting"
   - Quality: "8K detail", "sharp focus", "professional photography"

3. **Describe Composition:**
   - "centered composition", "rule of thirds", "symmetrical", "dynamic angle"
   - "foreground interest with blurred background"
   - "overhead flat lay composition"

4. **Specify Materials and Textures:**
   - "brushed metal texture", "soft fabric", "weathered wood grain"
   - "glossy finish", "matte surface", "translucent glass"

5. **Use Art References:**
   - "in the style of Studio Ghibli"
   - "impressionist painting technique"
   - "art nouveau design elements"

## Model Selection Guide

### Choose Flash When:
- Need rapid iteration (5-10+ variations)
- Working on concepts/drafts
- Output will be web/mobile resolution
- Budget is a concern
- Speed is priority over quality

### Choose Pro When:
- Creating final deliverables
- Need 2K or 4K resolution
- Image contains text that must be legible
- Composition is complex (10+ elements)
- Using multiple reference images
- Visualizing current events (with Search)
- Print quality required

## Cost Optimization

**Flash Model:**
- Use for all exploration and iteration
- Generate multiple concepts quickly
- Switch to Pro only for finals

**Pro Model:**
- Reserve for final selected concepts
- Use when resolution matters
- Leverage when complexity demands it

**Batch Processing:**
- Group similar requests
- Use chat context for refinements
- Save successful prompts for reuse

## Error Handling

### Common Issues

**Content Policy Violations:**
```python
try:
    response = client.models.generate_content(
        model="gemini-2.5-flash-image",
        contents=prompt,
        config=types.GenerateContentConfig(response_modalities=["TEXT", "IMAGE"])
    )
except Exception as e:
    if "SAFETY" in str(e):
        print("Content policy violation. Try rephrasing the prompt.")
    elif "QUOTA" in str(e):
        print("Rate limit reached. Wait before retrying.")
```

**Empty or Missing Images:**
```python
has_image = False
for part in response.parts:
    if part.inline_data is not None:
        has_image = True
        image = part.as_image()
        image.save("output.png")

if not has_image:
    print("No image generated. Try a more specific prompt.")
```

## Rate Limits and Quotas

- Free tier: Limited requests per minute/day
- Paid tier: Higher limits, consult current pricing
- Pro model: Lower rate limits than Flash
- Implement exponential backoff for retries

## Best Practices Summary

1. **Start with Flash** for exploration, move to Pro for finals
2. **Use chat interface** for iterative refinements
3. **Be specific** in prompts - detail generates better results
4. **Match aspect ratio** to intended use case
5. **Leverage Pro features** (Search, multiple refs) when needed
6. **Save intermediate results** - generation isn't deterministic
7. **Handle errors gracefully** - API limits and content policies apply
8. **Optimize costs** - use appropriate model for each task

## API Resources

- **SDK Migration Guide**: https://ai.google.dev/gemini-api/docs/migrate
- **Image Generation Docs**: https://ai.google.dev/gemini-api/docs/image-generation
- **API Reference**: https://ai.google.dev/api/generate-content
- **Python SDK (new)**: https://github.com/googleapis/python-genai
- **Deprecated SDK Info**: https://github.com/google-gemini/deprecated-generative-ai-python
- **Pricing**: https://ai.google.dev/pricing
- **Content Policies**: https://ai.google.dev/gemini-api/terms
